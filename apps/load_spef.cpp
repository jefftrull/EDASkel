// A simple SPEF loader
// Copyright (C) 2013 Jeffrey Elliot Trull <edaskel@att.net>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include <fstream>

// SPEF parser
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include "../parser/spefparser.hpp"

// Graph library for storing RC networks
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>

// vertex (circuit node) properties are just a string with the name we got out of SPEF;
typedef std::string vertex_property_t;

// edge (circuit component) properties are the type (resistor or capacitor) with the component value
using namespace boost::units;
typedef quantity<si::resistance, double> resistor_edge_t;
typedef quantity<si::capacitance, double> capacitor_edge_t;
typedef boost::variant<resistor_edge_t, capacitor_edge_t> edge_property_t;

typedef boost::adjacency_list<boost::vecS, boost::listS, boost::undirectedS,
                              vertex_property_t, edge_property_t> CktGraph;

struct Visitor {
  Visitor() : g(new CktGraph()), gnd(add_vertex(*g)) {}

  typedef size_t name_token_value_t;
  name_token_value_t name_map_entry(std::string n) {
    name_token_value_t id = names.size();
    names.push_back(std::move(n));
    return id;   // just the offset within the name map, for now
  }

  void port_definition(name_token_value_t /*net*/, char /*dir*/) {
  }

  void net_definition(name_token_value_t net,
                      quantity<si::capacitance, double> lumpc) {
    lumped_caps.emplace(net, lumpc);
  }

  void net_port_connection(name_token_value_t net, name_token_value_t port) {
    // add this to the list of connections for the net
    net_ports[net].push_back(port);   // creates map entry if not present
    get_vertex(port, "");    // ensure the vertex descriptor map contains this connection
  }

  void net_inst_connection(name_token_value_t net, name_token_value_t inst, std::string const& pin) {
    net_iconns[net].emplace_back(inst, pin);
    get_vertex(inst, pin);
  }

  void capacitor(name_token_value_t /*net*/, unsigned /*capnum*/,
                 name_token_value_t net_or_inst1, std::string const& node_or_pin1,
                 name_token_value_t net_or_inst2, std::string const& node_or_pin2,
                 quantity<si::capacitance, double> value) {
    // store in graph as edge
    boost::add_edge(get_vertex(net_or_inst1, node_or_pin1),
                    get_vertex(net_or_inst2, node_or_pin2),
                    value, *g);
  }

  void cgnd     (name_token_value_t /*net*/, unsigned /*capnum*/,
                 name_token_value_t net_or_inst, std::string const& node_or_pin,
                 quantity<si::capacitance, double> value) {
    boost::add_edge(get_vertex(net_or_inst, node_or_pin),
                    gnd,
                    value, *g);
  }

  void resistor (name_token_value_t /*net*/, unsigned /*capnum*/,
                 name_token_value_t net_or_inst1, std::string const& node_or_pin1,
                 name_token_value_t net_or_inst2, std::string const& node_or_pin2,
                 quantity<si::resistance, double> value) {
    boost::add_edge(get_vertex(net_or_inst1, node_or_pin1),
                    get_vertex(net_or_inst2, node_or_pin2),
                    value, *g);
  }

  std::vector<std::string> names;   // from the name map
  std::map<name_token_value_t, quantity<si::capacitance, double> > lumped_caps;

  std::shared_ptr<CktGraph> g;

  typedef std::map<std::pair<name_token_value_t, std::string>,
                   CktGraph::vertex_descriptor> VertexMap;
  VertexMap vertex_desc_map;

  typedef std::map<name_token_value_t, std::vector<name_token_value_t> > PortsMap;
  PortsMap net_ports;

  typedef std::vector<std::pair<name_token_value_t, std::string> > InstConnList;
  typedef std::map<name_token_value_t, InstConnList> InstConnsMap;
  InstConnsMap net_iconns;

  // will not allocate if not present (throws out_of_range, in fact)
  CktGraph::vertex_descriptor vertex(name_token_value_t net, std::string const& node_or_pin) const {
    return vertex_desc_map.at(std::make_pair(net, node_or_pin));
  }

private:
  CktGraph::vertex_descriptor gnd;

  // for internal use: allocates if not present
  CktGraph::vertex_descriptor get_vertex(name_token_value_t net, std::string const& node_or_pin) {
    // return index of net/node pair for use in circuit graph
    VertexMap::const_iterator it = vertex_desc_map.find(std::make_pair(net, node_or_pin));
    if (it == vertex_desc_map.end()) {
      bool found;  // it won't be
      std::tie(it, found) = vertex_desc_map.emplace(std::make_pair(net, node_or_pin),
                                                    add_vertex(*g)); // no vertex property for now
    }
    return it->second;
  }

};

// A Variant visitor for edge attributes that selects only resistors
struct IsResistor : boost::static_visitor<bool> {
  bool operator()(resistor_edge_t const&) const {
    return true;
  }
  bool operator()(capacitor_edge_t const&) const {
    return false;
  }
};

// a predicate class for creating a resistor-only filtered graph
struct ResistorsOnly {
  ResistorsOnly() {}
  ResistorsOnly(std::shared_ptr<CktGraph> g) : graph_(std::move(g)) {}
  bool operator()(CktGraph::edge_descriptor e) const {
    // access edge property (res/cap variant) and apply predicate visitor
    return boost::apply_visitor(IsResistor(), (*graph_)[e]);
  }
private:
  std::shared_ptr<CktGraph> graph_;
};

int main(int argc, char **argv) {
   using namespace std;
   if (argc < 1) {
      cerr << "usage: load_spef filename" << endl;
      return 1;
   }
   ifstream spefin(argv[1]);
   if (!spefin.is_open()) {
      cerr << "Failed to open SPEF file " << argv[1] << endl;
      return 1;
   }
   spefin.unsetf(ios::skipws);

   using namespace EDASkel::SpefParse;
   using boost::spirit::qi::phrase_parse;

   Visitor spefVisitor;
   spefparser<SpefIter, Visitor> spefParser(spefVisitor);
   spefskipper<SpefIter> spefSkipper;
   SpefIter beg(spefin), end;

   if (!phrase_parse(beg, end, spefParser, spefSkipper)) {
      cerr << "parsing failed" << endl;
      return 1;
   }
   if (beg != end) {
     cerr << "not all input processed.  Extra input: |";
     copy(beg, end, ostream_iterator<char>(cerr, ""));
     cerr << "|" << endl;
      return 1;
   }

   cout << "Parasitic graph contains " << num_vertices(*spefVisitor.g) << " vertices and " << num_edges(*spefVisitor.g) << " edges" << endl;

   // Produce a filtered graph containing only resistor edges
   ResistorsOnly res_filter(spefVisitor.g);
   typedef boost::filtered_graph<CktGraph, ResistorsOnly> ResGraph;
   ResGraph res_graph(*spefVisitor.g, res_filter);

   // Divide resistor-only circuit graph into connected components

   // define required color and component property maps
   typedef ResGraph::vertices_size_type comp_number_t;
   typedef map<ResGraph::vertex_descriptor, comp_number_t> CCompsStorageMap;
   CCompsStorageMap comps;      // component map for algorithm results
   typedef map<ResGraph::vertex_descriptor, boost::default_color_type> ColorsStorageMap;
   ColorsStorageMap colors;     // temp storage for algorithm
   
   // adapt std::map for use by Graph algorithms as "property map"
   boost::associative_property_map<CCompsStorageMap> cpmap(comps);
   boost::associative_property_map<ColorsStorageMap> clrpmap(colors);

   // Run the algorithm
   connected_components(res_graph, cpmap, boost::color_map(clrpmap));

   // At this point, "comps" contains, for each vertex, the component to which it belongs
   // We need to verify that the vertices associated with the connections (instance and port)
   // for each net are all in the same component

   // Iterate over the nets, verifying that all connections of each net are in the same connected component
   for (auto const& lcpair : spefVisitor.lumped_caps) {
      map<comp_number_t, vector<string> > comps2conn_names;  // every component and its attached connections

      // get component numbers for all port and instance connections
      auto net = lcpair.first;

      // Record components associated with ports
      auto npit = spefVisitor.net_ports.find(net);
      if (npit != spefVisitor.net_ports.end()) {
        // this net has port connections
        for (auto const& pname : npit->second) {
          // look up the vertex descriptor for this port
          ResGraph::vertex_descriptor port_desc = spefVisitor.vertex(pname, "");
          // add its component to our tracking info for this net
          comps2conn_names[comps.at(port_desc)].push_back(spefVisitor.names.at(pname));
        }
      }

      // components connected to instance pins
      auto niit = spefVisitor.net_iconns.find(net);
      if (niit != spefVisitor.net_iconns.end()) {
        // instance connections:
        for (auto const& ic : niit->second) {
          ResGraph::vertex_descriptor iconn_desc = spefVisitor.vertex(ic.first, ic.second);
          comps2conn_names[comps.at(iconn_desc)].push_back(spefVisitor.names.at(ic.first) + ":" + ic.second);
        }
      }

      if (comps2conn_names.size() > 1) {
        cout << "Net " << spefVisitor.names.at(net) << " appears to be discontinuous:" << endl;
        for (auto const& comp : comps2conn_names) {
          cout << "  Connection Group:" << endl;
          for (auto const& vtcomp : comp.second) {
             cout << "    " << vtcomp << endl;
          }
        }
      }
   }

}
