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

#include <iostream>
#include <fstream>

// SPEF parser
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include "../parser/spefparser.hpp"

using namespace boost::units;

struct Visitor {
  typedef size_t name_token_value_t;
  name_token_value_t name_map_entry(std::string n) {
    name_token_value_t id = names.size();
    names.push_back(n);
    return id;   // just the offset within the name map, for now
  }

  void port_definition(name_token_value_t net, char dir) {
  }

  void net_definition(name_token_value_t net,
                      quantity<si::capacitance, double> lumpc) {
    lumped_caps.emplace(net, lumpc);
  }

  void net_port_connection(name_token_value_t net, name_token_value_t port) {}

  void net_inst_connection(name_token_value_t net, name_token_value_t inst, std::string const& pin) {}

  void capacitor(name_token_value_t net, unsigned capnum,
                 name_token_value_t net_or_inst1, std::string const& node_or_pin1,
                 name_token_value_t net_or_inst2, std::string const& node_or_pin2,
                 quantity<si::capacitance, double> value) {}

  void cgnd     (name_token_value_t net, unsigned capnum,
                 name_token_value_t net_or_inst, std::string const& node_or_pin,
                 quantity<si::capacitance, double> value) {}

  void resistor (name_token_value_t net, unsigned capnum,
                 name_token_value_t net_or_inst1, std::string const& node_or_pin1,
                 name_token_value_t net_or_inst2, std::string const& node_or_pin2,
                 quantity<si::resistance, double> value) {}

  std::vector<std::string> names;
  std::map<name_token_value_t, quantity<si::capacitance, double> > lumped_caps;
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
     std::copy(beg, end, ostream_iterator<char>(cerr, ""));
     cerr << "|" << endl;
      return 1;
   }

   cout << "Total capacitance per net (from lumped cap information):" << endl;
   for (auto lc : spefVisitor.lumped_caps) {
      cout << lc.second << " " << spefVisitor.names.at(lc.first) << endl;
   }

}
