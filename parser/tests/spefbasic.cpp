// Basic tests for the SPEF parser
// Copyright (C) 2010 Jeffrey Elliot Trull <edaskel@att.net>
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

spefparser<SpefIter> spefParser;
spefskipper<SpefIter> spefSkipper;

// boilerplate parsing code
void parse_check(std::string const& str, spef& result) {
  std::stringstream testspef(str);
  testspef.unsetf(std::ios::skipws);
  SpefIter beg(testspef), end;
  BOOST_CHECK( phrase_parse(beg, end, spefParser, spefSkipper, result) );  // we should match
  BOOST_CHECK( (beg == end) );                        // we should consume all input
  if (beg != end) {
     std::cerr << "excess input: ";
     std::copy(beg, end, std::ostream_iterator<char>(std::cerr, ""));
     std::cerr << std::endl;
  }
}


void parse_check_fail(std::string const& str) {
  std::stringstream testspef(str);
  testspef.unsetf(std::ios::skipws);
  SpefIter beg(testspef), end;
  spef result;
  BOOST_CHECK( !phrase_parse(beg, end, spefParser, spefSkipper, result) );  // we should NOT match
}

BOOST_AUTO_TEST_CASE( design_name ) {

   spef result;
   std::string units("*T_UNIT 1 PS\n*C_UNIT 1 PF\n*R_UNIT 1 KOHM\n*L_UNIT 1 HENRY\n");
   std::string preamble("SPEF\n*SPEF \"IEEE 1481-1998\"\n*DESIGN \"GreatDesign\"\n");
   preamble += "*DIVIDER /\n*DELIMITER :\n*BUS_DELIMITER [ ]\n";

   parse_check(preamble + "// some other comment\n" + units, result);
   BOOST_CHECK_EQUAL( "GreatDesign", result.name );
   BOOST_CHECK_EQUAL( "IEEE 1481-1998", result.standard );

   parse_check_fail("// no initial SPEF\n*SPEF \"IEEE 1481-1998\"\n*DESIGN \"SomewhatLessGood\"\n" + units);
   parse_check_fail("SPEF\n// initial comment\n*SPEF \"IEEE 1481-1998\"\n// look, no design name here\n" + units);
   parse_check_fail("SPEF\n// some comment\n// missing standard string\n*DESIGN \"BadDesign\"\n// some other comment\n" + units);

}

BOOST_AUTO_TEST_CASE( units ) {

   spef result;
   std::string preamble("SPEF\n*SPEF \"IEEE 1481-1998\"\n*DESIGN \"GreatDesign\"\n");
   preamble += "*DIVIDER /\n*DELIMITER :\n*BUS_DELIMITER [ ]\n";
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
