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
  typedef size_t net_token_value_t;
  net_token_value_t name_map_entry(std::string n) {
    net_token_value_t id = names.size();
    names.push_back(n);
    return id;   // just the offset within the name map, for now
  }

  void port_definition(net_token_value_t net, char dir) {
    ports.emplace_back(net, dir);
  }

  void net_definition(net_token_value_t net, double lumpc) {
    lumped_caps.emplace(net, lumpc);
  }

  std::vector<std::string> names;
  typedef std::pair<net_token_value_t, char> port_t;
  std::vector<port_t> ports;
  std::map<net_token_value_t, double> lumped_caps;
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
  std::string preamble = spefData::preamble + spefData::empty_flow + spefData::delimiters;
  parse_check(preamble + "*T_UNIT 1 PS\n*C_UNIT 1 PF\n*R_UNIT 1 KOHM\n*L_UNIT 1 HENRY\n", result);
  BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-12 * si::seconds)), result.t_unit );
  BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1e3 * si::ohms)), result.r_unit );
  BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(1e-12 * si::farads)), result.c_unit );
  BOOST_CHECK_EQUAL( (quantity<si::inductance, double>(1 * si::henrys)), result.l_unit );

  parse_check(preamble + "*T_UNIT 1 NS\n*C_UNIT 1 FF\n*R_UNIT 1 KOHM\n*L_UNIT 1 UH\n", result);
  BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-9 * si::seconds)), result.t_unit );
  BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1e3 * si::ohms)), result.r_unit );
  BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(1e-15 * si::farads)), result.c_unit );
  BOOST_CHECK_EQUAL( (quantity<si::inductance, double>(1e-6 * si::henrys)), result.l_unit );

  parse_check(preamble + "*T_UNIT 1 US\n*C_UNIT 1 PF\n*R_UNIT 1 OHM\n*L_UNIT 1 HENRY\n", result);
  BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-6 * si::seconds)), result.t_unit );
  BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1 * si::ohms)), result.r_unit );

  parse_check(preamble + "*T_UNIT 1 MS\n*C_UNIT 2.0 PF\n*R_UNIT 1.0 KOHM\n*L_UNIT 1 HENRY\n", result);
  BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-3 * si::seconds)), result.t_unit );
  BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(2e-12 * si::farads)), result.c_unit );

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
  std::string nets_minimal("*D_NET *2 0.0011\n*END\n");

  Visitor spefVisitorMin;
  parse_check(spefData::header + name_map + nets_minimal, spefVisitorMin, result);
  BOOST_REQUIRE_EQUAL( 1, spefVisitorMin.lumped_caps.size() );
  BOOST_CHECK_CLOSE( 0.0011, spefVisitorMin.lumped_caps.at(1), 0.000001 );

  Visitor spefVisitor;
  std::string nets_unimplemented_features("*D_NET *1 29.33\n"
                                          "*CONN\n"
                                          "*CAP\n1 *2:2 99.21\n2 *3:1 0.88\n"
                                          "*RES\n1 *1:1 *1:2 9.8765\n"
                                          "*END\n");

  parse_check(spefData::header + name_map + nets_unimplemented_features,
              spefVisitor, result);
  BOOST_REQUIRE_EQUAL( 1, spefVisitor.lumped_caps.size() );
  BOOST_CHECK_CLOSE( 29.33, spefVisitor.lumped_caps.at(0), 0.000001 );

}
