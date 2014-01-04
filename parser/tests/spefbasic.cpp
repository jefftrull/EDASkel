// Basic tests for the SPEF parser
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

#define BOOST_TEST_MODULE basic SPEF tests
#include <boost/test/included/unit_test.hpp>

#include "../spefparser.hpp"
#include <string>
#include <iostream>
#include <sstream>

using namespace EDASkel::SpefParse;

using boost::spirit::qi::phrase_parse;

struct Visitor {
  typedef size_t name_token_value_t;
  name_token_value_t name_map_entry(std::string n) {
    name_token_value_t id = names.size();
    names.push_back(n);
    return id;   // just the offset within the name map, for now
  }

  void port_definition(name_token_value_t net, char dir) {
    ports.emplace_back(net, dir);
  }

  void net_definition(name_token_value_t net,
                      quantity<si::capacitance, double> lumpc) {
    total_caps.emplace(net, lumpc);
  }

  void net_port_connection(name_token_value_t net, name_token_value_t port) {
    if (port_connections.find(net) == port_connections.end()) {
      port_connections.insert(std::make_pair(net, std::vector<name_token_value_t>()));
    }
    port_connections[net].push_back(port);
  }

  void net_inst_connection(name_token_value_t net, name_token_value_t inst, std::string const& pin) {
    if (inst_connections.find(net) == inst_connections.end()) {
      inst_connections.insert(std::make_pair(net, std::vector<inst_connection_t>()));
    }
    inst_connections[net].push_back(std::make_pair(inst, pin));
  }

  void capacitor(name_token_value_t net, unsigned index,
                 name_token_value_t net_or_inst1, std::string const& node_or_pin1,
                 name_token_value_t net_or_inst2, std::string const& node_or_pin2,
                 quantity<si::capacitance, double> value) {
    using namespace std;

    // record this data for test assertions
    if (capacitors.find(net) == capacitors.end()) {
      capacitors.insert(make_pair(net, capacitance_map_t()));
    }
    BOOST_REQUIRE(capacitors.at(net).find(index) == capacitors.at(net).end());  // should be new
    capacitors.at(net).insert(make_pair(index,
                                        capacitor_t(make_pair(net_or_inst1, node_or_pin1),
                                                    make_pair(net_or_inst2, node_or_pin2),
                                                    value)));
  }

  void cgnd     (name_token_value_t net, unsigned index,
                 name_token_value_t net_or_inst, std::string const& node_or_pin,
                 quantity<si::capacitance, double> value) {
    using namespace std;

    if (gnd_lumped_caps.find(net) == gnd_lumped_caps.end()) {
      gnd_lumped_caps.insert(make_pair(net, capacitance_gnd_map_t()));
    }
    // Only one use of each parasitic index, per net
    BOOST_REQUIRE(gnd_lumped_caps.at(net).find(index) == gnd_lumped_caps.at(net).end());
    gnd_lumped_caps.at(net).insert(make_pair(index,
                                             capacitor_gnd_t(make_pair(net_or_inst, node_or_pin),
                                                             value)));
  }

  void resistor (name_token_value_t net, unsigned index,
                 name_token_value_t net_or_inst1, std::string const& node_or_pin1,
                 name_token_value_t net_or_inst2, std::string const& node_or_pin2,
                 quantity<si::resistance, double> value) {
    using namespace std;
    if (resistors.find(net) == resistors.end()) {
      resistors.insert(make_pair(net, resistance_map_t()));
    }
    BOOST_REQUIRE(resistors.at(net).find(index) == resistors.at(net).end());  // should be new
    resistors.at(net).insert(make_pair(index,
                                       resistor_t(make_pair(net_or_inst1, node_or_pin1),
                                                  make_pair(net_or_inst2, node_or_pin2),
                                                  value)));
  }

  // data stored in a manner convenient for testing.  Not efficient for apps...

  std::vector<std::string> names;
  typedef std::pair<name_token_value_t, char> port_t;
  std::vector<port_t> ports;
  std::map<name_token_value_t, quantity<si::capacitance, double> > total_caps;
  typedef std::pair<name_token_value_t, std::string> inst_connection_t;
  std::map<name_token_value_t, std::vector<inst_connection_t> > inst_connections;
  std::map<name_token_value_t, std::vector<name_token_value_t> > port_connections;

  typedef std::pair<name_token_value_t, std::string> connection_t;
  struct resistor_t {
    connection_t conn1;
    connection_t conn2;
    quantity<si::resistance, double> value;
    resistor_t(connection_t c1, connection_t c2, quantity<si::resistance, double> v)
      : conn1(c1), conn2(c2), value(v) {}
  };
  struct capacitor_t {
    connection_t conn1;
    connection_t conn2;
    quantity<si::capacitance, double> value;
    capacitor_t(connection_t c1, connection_t c2, quantity<si::capacitance, double> v)
      : conn1(c1), conn2(c2), value(v) {}
  };
  struct capacitor_gnd_t {
    connection_t conn;
    quantity<si::capacitance, double> value;
    capacitor_gnd_t(connection_t c, quantity<si::capacitance, double> v)
      : conn(c), value(v) {}
  };

  typedef std::map<unsigned, capacitor_gnd_t> capacitance_gnd_map_t;
  typedef std::map<unsigned, resistor_t>  resistance_map_t;
  typedef std::map<unsigned, capacitor_t> capacitance_map_t;
  typedef std::map<name_token_value_t, capacitance_gnd_map_t> net_lumpedcaps_t;
  typedef std::map<name_token_value_t, resistance_map_t> net_res_connections_t;
  typedef std::map<name_token_value_t, capacitance_map_t> cap_connections_t;
  net_lumpedcaps_t gnd_lumped_caps;
  net_res_connections_t resistors;
  cap_connections_t capacitors;

};

// without this we would have to make std::pair<size_t, char> iostream-able
BOOST_TEST_DONT_PRINT_LOG_VALUE(Visitor::port_t);

// boilerplate parsing code
void parse_check(std::string const& str, Visitor& spefVisitor, spef& result) {
  spefparser<SpefIter, Visitor> spefParser(spefVisitor);
  spefskipper<SpefIter> spefSkipper;

  std::stringstream testspef(str);
  testspef.unsetf(std::ios::skipws);
  SpefIter beg(testspef), end;
  BOOST_CHECK( phrase_parse(beg, end, spefParser, spefSkipper, result) );  // we should match
  BOOST_CHECK( (beg == end) );                        // we should consume all input

}

// for checks not needing data from the visitor
void parse_check(std::string const& str, spef& result) {
  Visitor spefVisitor;
  parse_check(str, spefVisitor, result);
}

void parse_check_fail(std::string const& str) {
  Visitor spefVisitor;
  spefparser<SpefIter, Visitor> spefParser(spefVisitor);
  spefskipper<SpefIter> spefSkipper;

  std::stringstream testspef(str);
  testspef.unsetf(std::ios::skipws);
  SpefIter beg(testspef), end;
  spef result;
  BOOST_CHECK( !phrase_parse(beg, end, spefParser, spefSkipper, result) );  // we should NOT match
}

namespace spefData {
  std::string preamble_no_magic("*SPEF \"IEEE 1481-1998\"\n*DESIGN \"GreatDesign\"\n"
                                "*DATE \"Fri Aug 29 02:14:00 1997\"\n*VENDOR \"MegaEDA\"\n"
                                "*PROGRAM \"Ozymandius Extractor\"\n*VERSION \"0.99\"\n");
  std::string preamble("SPEF\n" + preamble_no_magic);
  std::string empty_flow("*DESIGN_FLOW\n");
  std::string delimiters("*DIVIDER /\n*DELIMITER :\n*BUS_DELIMITER [ ]\n");
  std::string units("*T_UNIT 1 PS\n*C_UNIT 1 PF\n*R_UNIT 1 KOHM\n*L_UNIT 1 HENRY\n");
  std::string header(preamble + empty_flow + delimiters + units);

}

BOOST_AUTO_TEST_CASE( header ) {

  spef result;
  using namespace spefData;

  parse_check(preamble + "// some other comment\n" + empty_flow + delimiters + spefData::units, result);
  BOOST_CHECK_EQUAL( "GreatDesign", result.name );
  BOOST_CHECK_EQUAL( "IEEE 1481-1998", result.standard );
  BOOST_CHECK_EQUAL( "MegaEDA", result.vendor );
  BOOST_CHECK_EQUAL( "Ozymandius Extractor", result.program );
  BOOST_CHECK_EQUAL( "0.99", result.version );
  BOOST_CHECK(       result.design_flow.empty() );

  // verify we can omit the initial "SPEF"
  parse_check(preamble_no_magic + "// some other comment\n" + empty_flow + delimiters + spefData::units, result);

  std::string more_interesting_flow("*DESIGN_FLOW \"ALPHA BETA_GAMMA\" \"DELTA\"\n");
  parse_check(preamble + more_interesting_flow + delimiters + units, result);
  BOOST_CHECK_EQUAL( 2, result.design_flow.size() );
  BOOST_CHECK_EQUAL( "BETA_GAMMA", result.design_flow.find("ALPHA")->second );
  BOOST_CHECK_EQUAL( "", result.design_flow.find("DELTA")->second );

}

BOOST_AUTO_TEST_CASE( units ) {

  spef result;
  Visitor spefVisitor;

  std::string preamble = spefData::preamble + spefData::empty_flow + spefData::delimiters;
  std::string components = "*NAME_MAP\n*1 A\n*2 B\n"
                           "*D_NET *1 121.0\n"
                           "*CAP\n1 *1:1 1.0\n"
                           "*RES\n1 *1:1 *2:1 1.0\n"
                           "*END\n";

  parse_check(preamble + "*T_UNIT 1 PS\n*C_UNIT 1 PF\n*R_UNIT 1 KOHM\n*L_UNIT 1 HENRY\n" + components,
              spefVisitor, result);
  BOOST_REQUIRE_EQUAL(1, spefVisitor.resistors.size() );
  BOOST_REQUIRE_EQUAL(1, spefVisitor.gnd_lumped_caps.size() );
  BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-12 * si::seconds)), result.t_unit );
  BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1e3 * si::ohms)),
                     spefVisitor.resistors.at(0).at(1).value );
  BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(1e-12 * si::farads)),
                     spefVisitor.gnd_lumped_caps.at(0).at(1).value );
  BOOST_CHECK_EQUAL( (quantity<si::inductance, double>(1 * si::henrys)), result.l_unit );

  spefVisitor = Visitor();   // replace with fresh visitor
  parse_check(preamble + "*T_UNIT 1 NS\n*C_UNIT 1 FF\n*R_UNIT 1 KOHM\n*L_UNIT 1 UH\n" + components,
              spefVisitor, result);
  BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-9 * si::seconds)), result.t_unit );
  BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1e3 * si::ohms)),
                     spefVisitor.resistors.at(0).at(1).value );
  BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(1e-15 * si::farads)),
                     spefVisitor.gnd_lumped_caps.at(0).at(1).value );
  BOOST_CHECK_EQUAL( (quantity<si::inductance, double>(1e-6 * si::henrys)), result.l_unit );

  spefVisitor = Visitor();
  parse_check(preamble + "*T_UNIT 1 US\n*C_UNIT 1 PF\n*R_UNIT 1 OHM\n*L_UNIT 1 HENRY\n" + components,
              spefVisitor, result);
  BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-6 * si::seconds)), result.t_unit );
  BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1 * si::ohms)),
                     spefVisitor.resistors.at(0).at(1).value );

  spefVisitor = Visitor();
  parse_check(preamble + "*T_UNIT 1 MS\n*C_UNIT 2.0 PF\n*R_UNIT 1.0 KOHM\n*L_UNIT 1 HENRY\n" + components,
              spefVisitor, result);
  BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-3 * si::seconds)), result.t_unit );
  BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(2e-12 * si::farads)),
                     spefVisitor.gnd_lumped_caps.at(0).at(1).value );

  std::string crl("*C_UNIT 2.0 PF\n*R_UNIT 1.0 KOHM\n*L_UNIT 1 HENRY\n");
  parse_check_fail(preamble + "*T_UNIT 1 QS\n" + crl);  // a fake unit prefix
  parse_check_fail(preamble + "*X_UNIT 1 PS\n" + crl);  // no such thing as an "X" unit
  parse_check_fail(preamble + "*T_UNIT 1 PF\n" + crl);  // time but I've given capacitance

  // quantity and unit mismatches
  parse_check_fail(preamble + "*T_UNIT 1.0 PS\n*C_UNIT 1.0 PF\n*R_UNIT 1.0 FF\n*L_UNIT 1 HENRY\n");
  parse_check_fail(preamble + "*T_UNIT 1.0 PS\n*C_UNIT 1.0 OHM\n*R_UNIT 1.0 KOHM\n*L_UNIT 1 HENRY\n");

  // TODO: test for missing spaces? Not presently checked.
}

BOOST_AUTO_TEST_CASE( name_map ) {

  spef result;

  std::string no_map = spefData::header;
  parse_check(no_map, result);

  Visitor spefVisitor;
  parse_check(no_map + "*NAME_MAP\n*5 MOD1/BLK2/CellA[21]\n*10 top_lvl99\n",
              spefVisitor, result);

  BOOST_REQUIRE_EQUAL( 2, spefVisitor.names.size() );
  BOOST_CHECK_EQUAL( "MOD1/BLK2/CellA[21]", spefVisitor.names[0] );
  BOOST_CHECK_EQUAL( "top_lvl99", spefVisitor.names[1] );

  // check quoted delimiters
  Visitor spefVisitor2;
  parse_check(no_map + "*NAME_MAP\n*5 MOD1/BLK2/CellA\\[21\\]\n*10 top_lvl99\n",
              spefVisitor2, result);
  BOOST_REQUIRE_EQUAL( 2, spefVisitor2.names.size() );
  BOOST_CHECK_EQUAL( "MOD1/BLK2/CellA[21]", spefVisitor2.names[0] );
  BOOST_CHECK_EQUAL( "top_lvl99", spefVisitor2.names[1] );

}

BOOST_AUTO_TEST_CASE( ports ) {

  spef result;
  std::string name_map("*NAME_MAP\n*1 A/B/C\n*2 x1_12\n*3 z22[8]\n");

  std::string no_ports = spefData::header + name_map;
  Visitor spefVisitorNp;
  parse_check(no_ports, spefVisitorNp, result);
  BOOST_CHECK_EQUAL( 0, spefVisitorNp.ports.size() );

  // parse again, this time supplying some ports
  std::string with_ports(no_ports + "*PORTS\n" +
                         "*1 I *C 3.6 0\n" +
                         "*2 O *C 88.99 123.2\n" +
                         "*3 B *L 0.00 *S 0.0 0.0\n");
  Visitor spefVisitor;
  parse_check(with_ports, spefVisitor, result);

  BOOST_REQUIRE_EQUAL( 3, spefVisitor.ports.size() );
  BOOST_CHECK_EQUAL( (std::make_pair<size_t, char>(0, 'I')), spefVisitor.ports[0] );
  BOOST_CHECK_EQUAL( (std::make_pair<size_t, char>(1, 'O')), spefVisitor.ports[1] );
  BOOST_CHECK_EQUAL( (std::make_pair<size_t, char>(2, 'B')), spefVisitor.ports[2] );
   
}

BOOST_AUTO_TEST_CASE( nets ) {

  spef result;
  std::string name_map("*NAME_MAP\n*1 A/B/C\n*2 x1_12\n*3 z22[8]\n");
  std::string nets_minimal("*D_NET *2 0.0011\n*CONN\n*P *2 I\n*I *1:S I\n*END\n");

  Visitor spefVisitorMin;
  parse_check(spefData::header + name_map + nets_minimal, spefVisitorMin, result);
  BOOST_REQUIRE_EQUAL( 1, spefVisitorMin.total_caps.size() );
  BOOST_CHECK_CLOSE( 0.0011e-12, spefVisitorMin.total_caps.at(1).value(), 0.00001 );
  BOOST_REQUIRE_EQUAL( 1, spefVisitorMin.inst_connections.size() );
  // key of first entry is the name of the net
  Visitor::name_token_value_t firstNet = spefVisitorMin.inst_connections.begin()->first;
  BOOST_CHECK_EQUAL( "x1_12", spefVisitorMin.names[firstNet] );
  // get first entry of vector associated with the first key.
  // The first value of that pair is the name of the instance
  Visitor::name_token_value_t firstInst = spefVisitorMin.inst_connections.begin()->second.begin()->first;
  BOOST_CHECK_EQUAL( "A/B/C", spefVisitorMin.names[firstInst] );
  // The second value of that pair is the name of the pin
  BOOST_CHECK_EQUAL( "S", spefVisitorMin.inst_connections.begin()->second.begin()->second);

  Visitor spefVisitor;
  std::string nets_unimplemented_features("*D_NET *1 29.33\n"
                                          "*CAP\n1 *2:2 99.21\n2 *3:1 0.88\n"
                                          "*RES\n1 *1:1 *1:2 9.8765\n"
                                          "*END\n");

  parse_check(spefData::header + name_map + nets_unimplemented_features,
              spefVisitor, result);
  BOOST_REQUIRE_EQUAL( 1, spefVisitor.total_caps.size() );
  BOOST_CHECK_CLOSE( 29.33e-12, spefVisitor.total_caps.at(0).value(), 0.00001 );

}

BOOST_AUTO_TEST_CASE( simple_circuit ) {
  // construct a simple circuit consisting of two inputs connected to a nand,
  // followed by an inverter connected to an output, with parasitics in between
  std::string name_map("*NAME_MAP\n"
                       "*1 I1\n"
                       "*2 I2\n"
                       "*3 O_X\n"   // nand2 output
                       "*4 O\n"     // circuit output port
                       "*5 U1\n"    // nand2
                       "*6 U2\n");  // output inverter

  std::string ports("*PORTS\n"
                    "*1 I\n"
                    "*2 I\n"
                    "*4 O\n");

  std::string i1_net("*D_NET *1 0\n" // fake lumped C
                     "*CONN\n"
                     "*P *1 I\n"
                     "*I *5:A I\n"
                     "*CAP\n"            // 2 pi model connects input to pin
                     "1 *1:1 0.1\n"      // first cap of first pi model
                     "2 *1:2 0.2\n"      // middle two pi model caps
                     "3 *1:2 *2:2 0.1\n" // coupling cap between inputs
                     "4 *5:A 0.1\n"      // final pi model final cap
                     "*RES\n"
                     "1 *1:1 *1:2 100\n"
                     "2 *1:2 *5:A 100\n"
                     "*END\n");

  std::string i2_net("*D_NET *2 0\n"     // fake lumped C
                     "*CONN\n"
                     "*P *2 I\n"
                     "*I *5:B I\n"
                     "*CAP\n"            // identical to I1 net
                     "1 *2:1 0.1\n"
                     "2 *2:2 0.2\n"
                     "3 *2:2 *1:2 0.1\n"
                     "4 *5:B 0.1\n"
                     "*RES\n"
                     "1 *2:1 *2:2 100\n"
                     "2 *2:2 *5:B 100\n"
                     "*END\n");

  std::string ox_net("*D_NET *3 0.2\n"   // lumped C reflecting wire only(?)
                     "*CONN\n"
                     "*I *5:Z O\n"
                     "*I *6:A I\n"
                     "*CAP\n"
                     "1 *5:Z 0.1\n"      // a single pi model
                     "2 *6:A 0.1\n"
                     "*RES\n"
                     "1 *5:Z *6:A 100\n"
                     "*END\n");

  std::string o_net("*D_NET *4 0\n"
                     "*CONN\n"
                     "*P *4 O\n"
                     "*I *6:Z O\n"
                     "*CAP\n"          
                     "1 *6:Z 0.1\n"      // a simple lumped C
                     "*END\n");

  spef result;
  Visitor spefVisitor;
  parse_check(spefData::header + name_map + i1_net + i2_net + ox_net + o_net,
              spefVisitor, result);

  // check coupling cap
  BOOST_CHECK_EQUAL(2, spefVisitor.capacitors.size());  // only one net-to-net cap; it appears twice
  // I1/I2 are the first two entries in the name map and so should be indexed that way:
  BOOST_REQUIRE(spefVisitor.capacitors.find(0) != spefVisitor.capacitors.end());
  BOOST_REQUIRE(spefVisitor.capacitors.find(1) != spefVisitor.capacitors.end());
  // parasitic index for both is 3:
  BOOST_REQUIRE(spefVisitor.capacitors.at(0).find(3) != spefVisitor.capacitors.at(0).end());
  BOOST_REQUIRE(spefVisitor.capacitors.at(1).find(3) != spefVisitor.capacitors.at(1).end());
  // check connections
  // I1 first:
  auto i1_conn = spefVisitor.capacitors.at(0).at(3);
  BOOST_CHECK_EQUAL(0, i1_conn.conn1.first);         // "me" net
  BOOST_CHECK_EQUAL("2", i1_conn.conn1.second);      // my node
  BOOST_CHECK_EQUAL(1, i1_conn.conn2.first);         // other net
  BOOST_CHECK_EQUAL("2", i1_conn.conn2.second);      // other node
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), i1_conn.value);
  // then I2
  auto i2_conn = spefVisitor.capacitors.at(1).at(3);
  BOOST_CHECK_EQUAL(1, i2_conn.conn1.first);
  BOOST_CHECK_EQUAL("2", i2_conn.conn1.second);
  BOOST_CHECK_EQUAL(0, i2_conn.conn2.first);
  BOOST_CHECK_EQUAL("2", i2_conn.conn2.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), i2_conn.value);

  // check grounded caps
  BOOST_CHECK_EQUAL(4, spefVisitor.gnd_lumped_caps.size());  // 4 nets have cap to gnd

  // I1 has 3 (beginning, intermediate (x2), and final, from pi models
  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.find(0) != spefVisitor.gnd_lumped_caps.end());
  BOOST_CHECK_EQUAL(3, spefVisitor.gnd_lumped_caps.at(0).size());

  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(0).find(1) != spefVisitor.gnd_lumped_caps.at(0).end());
  auto i1_c1 = spefVisitor.gnd_lumped_caps.at(0).at(1);
  BOOST_CHECK_EQUAL(0, i1_c1.conn.first);
  BOOST_CHECK_EQUAL("1", i1_c1.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), i1_c1.value);

  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(0).find(2) != spefVisitor.gnd_lumped_caps.at(0).end());
  auto i1_c2 = spefVisitor.gnd_lumped_caps.at(0).at(2);
  BOOST_CHECK_EQUAL(0, i1_c2.conn.first);
  BOOST_CHECK_EQUAL("2", i1_c2.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.2e-12 * si::farads)), i1_c2.value);

  // coupling cap is index 3, not present here
  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(0).find(4) != spefVisitor.gnd_lumped_caps.at(0).end());
  auto i1_c3 = spefVisitor.gnd_lumped_caps.at(0).at(4);
  BOOST_CHECK_EQUAL(4, i1_c3.conn.first);      // A pin of nand2
  BOOST_CHECK_EQUAL("A", i1_c3.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), i1_c3.value);

  // I2 is the same as I1
  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.find(1) != spefVisitor.gnd_lumped_caps.end());
  BOOST_CHECK_EQUAL(3, spefVisitor.gnd_lumped_caps.at(1).size());

  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(1).find(1) != spefVisitor.gnd_lumped_caps.at(1).end());
  auto i2_c1 = spefVisitor.gnd_lumped_caps.at(1).at(1);
  BOOST_CHECK_EQUAL(1, i2_c1.conn.first);
  BOOST_CHECK_EQUAL("1", i2_c1.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), i2_c1.value);

  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(1).find(2) != spefVisitor.gnd_lumped_caps.at(1).end());
  auto i2_c2 = spefVisitor.gnd_lumped_caps.at(1).at(2);
  BOOST_CHECK_EQUAL(1, i2_c2.conn.first);
  BOOST_CHECK_EQUAL("2", i2_c2.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.2e-12 * si::farads)), i2_c2.value);

  // coupling cap is index 3, not present here
  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(1).find(4) != spefVisitor.gnd_lumped_caps.at(1).end());
  auto i2_c3 = spefVisitor.gnd_lumped_caps.at(1).at(4);
  BOOST_CHECK_EQUAL(4, i2_c3.conn.first);      // B pin of nand2
  BOOST_CHECK_EQUAL("B", i2_c3.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), i2_c3.value);

  // O_X (nand2 output/inverter input)
  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.find(2) != spefVisitor.gnd_lumped_caps.end());
  BOOST_CHECK_EQUAL(2, spefVisitor.gnd_lumped_caps.at(2).size());  // single pi model

  // first lumped cap (first node of pi model)
  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(2).find(1) != spefVisitor.gnd_lumped_caps.at(2).end());
  auto ox_c1 = spefVisitor.gnd_lumped_caps.at(2).at(1);
  BOOST_CHECK_EQUAL(4, ox_c1.conn.first);      // Z pin of nand2
  BOOST_CHECK_EQUAL("Z", ox_c1.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), ox_c1.value);

  // second lumped cap (last node of pi model)
  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(2).find(2) != spefVisitor.gnd_lumped_caps.at(2).end());
  auto ox_c2 = spefVisitor.gnd_lumped_caps.at(2).at(2);
  BOOST_CHECK_EQUAL(5, ox_c2.conn.first);      // A pin of inverter
  BOOST_CHECK_EQUAL("A", ox_c2.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), ox_c2.value);

  // lumped cap at inverter output
  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.find(3) != spefVisitor.gnd_lumped_caps.end());
  BOOST_CHECK_EQUAL(1, spefVisitor.gnd_lumped_caps.at(3).size());

  BOOST_REQUIRE(spefVisitor.gnd_lumped_caps.at(3).find(1) != spefVisitor.gnd_lumped_caps.at(3).end());
  auto o_c = spefVisitor.gnd_lumped_caps.at(3).at(1);
  BOOST_CHECK_EQUAL(5, o_c.conn.first);      // Z pin of inverter
  BOOST_CHECK_EQUAL("Z", o_c.conn.second);
  BOOST_CHECK_EQUAL((quantity<si::capacitance, double>(0.1e-12 * si::farads)), o_c.value);

  // Now the resistors
  // They are present in the first three nets but not the last
  BOOST_CHECK_EQUAL(3, spefVisitor.resistors.size());

  // I1
  BOOST_REQUIRE(spefVisitor.resistors.find(0) != spefVisitor.resistors.end());
  BOOST_REQUIRE(spefVisitor.resistors.at(0).find(1) != spefVisitor.resistors.at(0).end());
  auto i1_r1 = spefVisitor.resistors.at(0).at(1);
  BOOST_CHECK_EQUAL(0, i1_r1.conn1.first);    // net I1
  BOOST_CHECK_EQUAL("1", i1_r1.conn1.second); // starts at node 1
  BOOST_CHECK_EQUAL(0, i1_r1.conn2.first);    // net I1
  BOOST_CHECK_EQUAL("2", i1_r1.conn2.second); // ends at node 2
  BOOST_CHECK_EQUAL((quantity<si::resistance, double>(100e3 * si::ohms)),
                    i1_r1.value);             // 100 ohms
  BOOST_REQUIRE(spefVisitor.resistors.at(0).find(2) != spefVisitor.resistors.at(0).end());
  auto i1_r2 = spefVisitor.resistors.at(0).at(2);
  BOOST_CHECK_EQUAL(0, i1_r2.conn1.first);    // net I1
  BOOST_CHECK_EQUAL("2", i1_r2.conn1.second); // starts at node 2
  BOOST_CHECK_EQUAL(4, i1_r2.conn2.first);    // instance U1
  BOOST_CHECK_EQUAL("A", i1_r2.conn2.second); // ends at pin A
  BOOST_CHECK_EQUAL((quantity<si::resistance, double>(100e3 * si::ohms)),
                    i1_r2.value);        // 100 ohms

  // I2
  BOOST_REQUIRE(spefVisitor.resistors.find(1) != spefVisitor.resistors.end());
  BOOST_REQUIRE(spefVisitor.resistors.at(1).find(1) != spefVisitor.resistors.at(1).end());
  auto i2_r1 = spefVisitor.resistors.at(1).at(1);
  BOOST_CHECK_EQUAL(1, i2_r1.conn1.first);    // net I2
  BOOST_CHECK_EQUAL("1", i2_r1.conn1.second); // starts at node 1
  BOOST_CHECK_EQUAL(1, i2_r1.conn2.first);    // net I2
  BOOST_CHECK_EQUAL("2", i2_r1.conn2.second); // ends at node 2
  BOOST_CHECK_EQUAL((quantity<si::resistance, double>(100e3 * si::ohms)),
                    i2_r1.value);        // 100 ohms
  BOOST_REQUIRE(spefVisitor.resistors.at(1).find(2) != spefVisitor.resistors.at(1).end());
  auto i2_r2 = spefVisitor.resistors.at(1).at(2);
  BOOST_CHECK_EQUAL(1, i2_r2.conn1.first);    // net I2
  BOOST_CHECK_EQUAL("2", i2_r2.conn1.second); // starts at node 2
  BOOST_CHECK_EQUAL(4, i2_r2.conn2.first);    // instance U1
  BOOST_CHECK_EQUAL("B", i2_r2.conn2.second); // ends at pin B
  BOOST_CHECK_EQUAL((quantity<si::resistance, double>(100e3 * si::ohms)),
                    i2_r2.value);        // 100 ohms

  // O_X
  BOOST_REQUIRE(spefVisitor.resistors.find(1) != spefVisitor.resistors.end());
  BOOST_REQUIRE(spefVisitor.resistors.at(2).find(1) != spefVisitor.resistors.at(2).end());
  auto ox_r = spefVisitor.resistors.at(2).at(1);
  BOOST_CHECK_EQUAL(4, ox_r.conn1.first);     // U1 pin Z
  BOOST_CHECK_EQUAL("Z", ox_r.conn1.second);
  BOOST_CHECK_EQUAL(5, ox_r.conn2.first);     // to U2 pin A
  BOOST_CHECK_EQUAL("A", ox_r.conn2.second);
  BOOST_CHECK_EQUAL((quantity<si::resistance, double>(100e3 * si::ohms)),
                    ox_r.value);

}
