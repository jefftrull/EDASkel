// Basic tests for the SPEF parser

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
   std::string units("*T_UNIT 1 PS\n*R_UNIT 1KOHM\n*C_UNIT 1PF\n");
   std::string preamble("SPEF\n*SPEF \"IEEE 1481-1998\"\n*DESIGN \"GreatDesign\"\n");
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
   parse_check(preamble + "*T_UNIT 1.0 PS\n*R_UNIT 1.0 KOHM\n*C_UNIT 1.0 PF\n", result);
   BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-12 * si::seconds)), result.t_unit );
   BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1e3 * si::ohms)), result.r_unit );
   BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(1e-12 * si::farads)), result.c_unit );

   parse_check(preamble + "*T_UNIT 1.0 NS\n*R_UNIT 1.0 KOHM\n*C_UNIT 1.0 FF\n", result);
   BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-9 * si::seconds)), result.t_unit );
   BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1e3 * si::ohms)), result.r_unit );
   BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(1e-15 * si::farads)), result.c_unit );

   parse_check(preamble + "*T_UNIT 1.0 US\n*R_UNIT 1.0 OHM\n*C_UNIT 1.0 PF\n", result);
   BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-6 * si::seconds)), result.t_unit );
   BOOST_CHECK_EQUAL( (quantity<si::resistance, double>(1 * si::ohms)), result.r_unit );

   parse_check(preamble + "*T_UNIT 1.0 MS\n*R_UNIT 1.0 KOHM\n*C_UNIT 2.0 PF\n", result);
   BOOST_CHECK_EQUAL( (quantity<si::time, double>(1e-3 * si::seconds)), result.t_unit );
   BOOST_CHECK_EQUAL( (quantity<si::capacitance, double>(2e-12 * si::farads)), result.c_unit );

   parse_check_fail(preamble + "*T_UNIT 1 QS\n");  // a fake unit prefix
   parse_check_fail(preamble + "*X_UNIT 1 PS\n");  // no such thing as an "X" unit
   parse_check_fail(preamble + "*T_UNIT 1 PF\n");  // time but I've given capacitance

   // quantity and unit mismatches
   parse_check_fail(preamble + "*T_UNIT 1.0 PS\n*R_UNIT 1.0 FF\n*C_UNIT 1.0 PF\n");
   parse_check_fail(preamble + "*T_UNIT 1.0 PS\n*R_UNIT 1.0 KOHM\n*C_UNIT 1.0 OHM\n");
}
