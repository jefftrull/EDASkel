// A Boost Spirit-based DEF parser, part of EDASkel, a sample EDA app
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

// This file contains some simple unit tests for the DEF parser
#define BOOST_TEST_MODULE basic test
#include <boost/test/included/unit_test.hpp>

#include "../defparser.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace DefParse;
namespace EDASkel {
  extern defparser<LefDefIter> defParser;
  extern lefdefskipper<LefDefIter> lefdefSkipper;

}
using namespace EDASkel;
using boost::spirit::qi::phrase_parse;

// boilerplate parsing code
void parse_check(std::string const& str, def& result) {
  std::stringstream testdef(str);
  testdef.unsetf(std::ios::skipws);
  LefDefIter beg(testdef), end;
  BOOST_CHECK( phrase_parse(beg, end, defParser, lefdefSkipper, result) );  // we should match
  BOOST_CHECK( (beg == end) );                        // we should consume all input
}

void parse_check_fail(std::string const& str) {
  std::stringstream testdef(str);
  testdef.unsetf(std::ios::skipws);
  LefDefIter beg(testdef), end;
  def result;
  BOOST_CHECK( !phrase_parse(beg, end, defParser, lefdefSkipper, result) );  // we should NOT match
}

BOOST_AUTO_TEST_CASE( version_parse_simple ) {

  def result;
  parse_check("DESIGN test ;\nVERSION 1.211 ;\nEND DESIGN\n", result);

  BOOST_CHECK_EQUAL( result.name, "test" );
  BOOST_CHECK_EQUAL( result.version, 1.211 );
}

BOOST_AUTO_TEST_CASE ( version_parse_nospace ) {
  // should fail ("distinct" issue)
  parse_check_fail("DESIGN test ;\nVERSION1.211 ;\nEND DESIGN\n");

}

BOOST_AUTO_TEST_CASE ( version_parse_nonnum ) {
  parse_check_fail("DESIGN test ;\nVERSION 1.21a ;\nEND DESIGN\n");

}

BOOST_AUTO_TEST_CASE ( version_parse_spaced_keywd ) {
  parse_check_fail("DESIGN test ;\nVER SION 1.211 ;\nEND DESIGN\n");

}

BOOST_AUTO_TEST_CASE ( components_parse_empty ) {
  def result;
  parse_check("DESIGN test ;\nCOMPONENTS 0 ;\nEND COMPONENTS\nEND DESIGN\n", result);

  BOOST_CHECK( result.components.empty() );
}
  
BOOST_AUTO_TEST_CASE ( components_parse_simple ) {
  def result;
  parse_check("DESIGN test-hyphenated ;\nVERSION 1.211 ;\nDIEAREA ( 0 0 ) ( 100000 200000 ) ;\nCOMPONENTS 1 ;\n - I111_uscore/hiername INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nSITE CORE1 0 0 N DO 200 BY 1 STEP 100 500 ;\nEND DESIGN\n", result);

  BOOST_CHECK_EQUAL( result.name, "test-hyphenated" );
  BOOST_CHECK_EQUAL( result.diearea.ll.x, 0 );
  BOOST_CHECK_EQUAL( result.diearea.ll.y, 0 );
  BOOST_CHECK_EQUAL( result.diearea.ur.x, 100000 );
  BOOST_CHECK_EQUAL( result.diearea.ur.y, 200000 );
  BOOST_REQUIRE_EQUAL( result.components.size(), 1 );    // exactly one component read
  BOOST_CHECK_EQUAL( result.components[0].name, "I111_uscore/hiername" );
  BOOST_CHECK_EQUAL( result.components[0].celltype, "INVX2" );
  BOOST_REQUIRE( (result.components[0].placement) );     // it has placement
  BOOST_CHECK_EQUAL( result.components[0].placement->plcfix, "FIXED" );
  BOOST_CHECK_EQUAL( result.components[0].placement->origin.x, -4107 );
  BOOST_CHECK_EQUAL( result.components[0].placement->origin.y, 82000 );
  BOOST_CHECK_EQUAL( result.components[0].placement->orient, "FN" );
}
  
BOOST_AUTO_TEST_CASE ( components_noplace ) {
  def result;
  parse_check("DESIGN test ;\nCOMPONENTS 1 ;\n - I111 INVX2 ;\nEND COMPONENTS\nEND DESIGN\n", result);

  BOOST_CHECK_EQUAL( result.name, "test" );
  BOOST_REQUIRE_EQUAL( result.components.size(), 1 );    // exactly one component read
  BOOST_CHECK_EQUAL( result.components[0].name, "I111" );
  BOOST_CHECK_EQUAL( result.components[0].celltype, "INVX2" );
  BOOST_CHECK( (!result.components[0].placement) );     // it has NO placement
}  

BOOST_AUTO_TEST_CASE ( components_parse_wrongcount ) {
  parse_check_fail("DESIGN test ;\nCOMPONENTS 2 ;\n - I111 INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nEND DESIGN\n");

}
  
BOOST_AUTO_TEST_CASE ( site_basic ) {
  def result;
  parse_check("DESIGN test ;\nVERSION 1.211 ;\nDIEAREA ( 0 0 ) ( 100000 200000 ) ;\nSITE CORE1 10 20 N DO 200 BY 1 STEP 100 500 ;\nEND DESIGN\n", result);

  BOOST_CHECK_EQUAL( result.diearea.ll.x, 0 );
  BOOST_CHECK_EQUAL( result.diearea.ll.y, 0 );
  BOOST_CHECK_EQUAL( result.diearea.ur.x, 100000 );
  BOOST_CHECK_EQUAL( result.diearea.ur.y, 200000 );
  BOOST_REQUIRE_EQUAL( result.rows.size(), 1 );
  BOOST_CHECK( !result.rows[0].rowname );
  BOOST_CHECK_EQUAL( result.rows[0].body.sitename, "CORE1" );
  BOOST_CHECK_EQUAL( result.rows[0].body.origx, 10 );
  BOOST_CHECK_EQUAL( result.rows[0].body.origy, 20 );
  BOOST_CHECK_EQUAL( result.rows[0].body.orient, "N" );
  BOOST_REQUIRE( result.rows[0].body.repeat );
  BOOST_CHECK_EQUAL( result.rows[0].body.repeat->xrepeat, 200 );
  BOOST_CHECK_EQUAL( result.rows[0].body.repeat->yrepeat, 1 );
  BOOST_REQUIRE( result.rows[0].body.repeat->step );
  BOOST_CHECK_EQUAL( result.rows[0].body.repeat->step->first, 100 );
  BOOST_CHECK_EQUAL( result.rows[0].body.repeat->step->second, 500 );
}

// BOZO when we eventually parse everything this won't be a very interesting test and probably should be removed,
// or have previously ignored stuff checked
BOOST_AUTO_TEST_CASE ( parse_ignored_stuff ) {
  def result;
  parse_check("DESIGN test ;\nVERSION 1.211 ;\nDIEAREA ( 0 0 ) ( 100000 200000 ) ;\nCOMPONENTS 1 ;\n - I111 INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nSITE CORE1 0 0 N DO 200 BY 1 STEP 100 500 ;\nSPECIALNETS 1 ;\n - GND ;\nEND SPECIALNETS\nEND DESIGN\n", result);

  BOOST_CHECK_EQUAL( result.name, "test" );
  BOOST_CHECK_EQUAL( result.diearea.ll.x, 0 );
  BOOST_CHECK_EQUAL( result.diearea.ll.y, 0 );
  BOOST_CHECK_EQUAL( result.diearea.ur.x, 100000 );
  BOOST_CHECK_EQUAL( result.diearea.ur.y, 200000 );
  BOOST_REQUIRE_EQUAL( result.components.size(), 1 );    // exactly one component read
  BOOST_CHECK_EQUAL( result.components[0].name, "I111" );
  BOOST_CHECK_EQUAL( result.components[0].celltype, "INVX2" );
  BOOST_REQUIRE( (result.components[0].placement) );     // it has placement
  BOOST_CHECK_EQUAL( result.components[0].placement->plcfix, "FIXED" );
  BOOST_CHECK_EQUAL( result.components[0].placement->origin.x, -4107 );
  BOOST_CHECK_EQUAL( result.components[0].placement->origin.y, 82000 );
  BOOST_CHECK_EQUAL( result.components[0].placement->orient, "FN" );
}
  
BOOST_AUTO_TEST_CASE( history ) {
  // history entries can have any random thing, some of which resembles other valid syntax
  def result;
  parse_check("DESIGN test ;\nHISTORY basic history line;\nHISTORY extra keywords VERSION 1.211 DIEAREA ( 0 0 ) ( 10 10 ) COMPONENTS END COMPONENTS;\nHISTORY SPECIALNETS 1 more keywords END SPECIALNETS;\nEND DESIGN\n", result);

  BOOST_CHECK_EQUAL( result.name, "test" );
  BOOST_REQUIRE_EQUAL( 3, result.history.size() );
  BOOST_CHECK_EQUAL( "basic history line", result.history[0] );
  BOOST_CHECK_EQUAL( "extra keywords VERSION 1.211 DIEAREA ( 0 0 ) ( 10 10 ) COMPONENTS END COMPONENTS", result.history[1] );
  BOOST_CHECK_EQUAL( "SPECIALNETS 1 more keywords END SPECIALNETS", result.history[2] );
}

BOOST_AUTO_TEST_CASE( net_simple ) {
  def result;
  parse_check("DESIGN test ;\nCOMPONENTS 2;\n- X C ;\n- Y C ;\nEND COMPONENTS\nNETS 3 ;\n- ALPHA ( X P1 ) ( Y P2 ) ;\n- BETA ;\n- GAMMA ;\nEND NETS\nEND DESIGN\n", result);

  BOOST_CHECK_EQUAL( result.name, "test" );
  BOOST_REQUIRE_EQUAL( 3, result.nets.size() );
  BOOST_CHECK_EQUAL( "ALPHA", result.nets[0].name );
  BOOST_CHECK_EQUAL( "BETA", result.nets[1].name );
  BOOST_CHECK_EQUAL( "GAMMA", result.nets[2].name );
  BOOST_REQUIRE_EQUAL( 2, result.nets[0].connections.size() );
  BOOST_CHECK_EQUAL( "X", result.nets[0].connections[0].first );
  BOOST_CHECK_EQUAL( "P1", result.nets[0].connections[0].second );
  BOOST_CHECK_EQUAL( "Y", result.nets[0].connections[1].first );
  BOOST_CHECK_EQUAL( "P2", result.nets[0].connections[1].second );
}

BOOST_AUTO_TEST_CASE( net_wrong_count ) {
  // both too many and too few
  parse_check_fail("DESIGN test ;\nNETS 2 ;\n- ALPHA ;\n- BETA ;\n- GAMMA ;\nEND NETS\nEND DESIGN\n");
  parse_check_fail("DESIGN test ;\nNETS 2 ;\n- ALPHA ;\nEND NETS\nEND DESIGN\n");
}

BOOST_AUTO_TEST_CASE( erroneous_connection ) {
  // no components, but the one net refers to a pin on one
  parse_check_fail("DESIGN test ;\nCOMPONENTS 0;\nEND COMPONENTS\nNETS 1 ;\n- OMEGA ( C1 A ) ;\nEND NETS\nEND DESIGN\n");
  // two components, but one net refers to a nonexistent component
  parse_check_fail("DESIGN test ;\nCOMPONENTS 2;\n- X C ;\n- Y C ;\nEND COMPONENTS\nNETS 3 ;\n- ALPHA ( BOGUS P1 ) ( Y P2 ) ;\n- BETA ( Y SOMEPIN ) ;\n- GAMMA ;\nEND NETS\nEND DESIGN\n");
}
