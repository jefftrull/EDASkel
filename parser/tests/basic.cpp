// A Boost Spirit-based DEF parser, part of EDASkel, a sample EDA app
// Copyright (C) 2010 Jeffrey Elliot Trull <linmodemstudent@gmail.com>
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

using namespace DefParse;
defparser<std::string::const_iterator> defp;
using boost::spirit::qi::space;

BOOST_AUTO_TEST_CASE( version_parse_simple ) {

  std::string testdef("DESIGN test ;\nVERSION 1.211 ;\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  BOOST_CHECK( phrase_parse(beg, end, defp, space) );  // we should match
  BOOST_CHECK( (beg == end) );                         // we should consume all input
  // did not use BOOST_CHECK_EQUAL b/c it wants to output these on failure, and there is no operator<< defined

}

BOOST_AUTO_TEST_CASE ( version_parse_nospace ) {

  std::string testdef("DESIGN test ;\nVERSION1.211 ;\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  BOOST_CHECK( !phrase_parse(beg, end, defp, space) );

}

BOOST_AUTO_TEST_CASE ( version_parse_nonnum ) {

  std::string testdef("DESIGN test ;\nVERSION 1.21a ;\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  BOOST_CHECK( !phrase_parse(beg, end, defp, space) );

}

BOOST_AUTO_TEST_CASE ( components_parse_empty ) {
  std::string testdef("DESIGN test ;\nCOMPONENTS 0 ;\nEND COMPONENTS\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  def result;
  BOOST_CHECK( phrase_parse(beg, end, defp, space, result) );
  BOOST_CHECK( (beg == end) );                         // we should consume all input

  BOOST_CHECK( result.components.empty() );
}
  
BOOST_AUTO_TEST_CASE ( components_parse_simple ) {
  std::string testdef("DESIGN test-hyphenated ;\nVERSION 1.211 ;\nDIEAREA ( 0 0 ) ( 100000 200000 ) ;\nCOMPONENTS 1 ;\n - I111_uscore/hiername INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nSITE CORE1 0 0 N DO 200 BY 1 STEP 100 500 ;\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  def result;
  BOOST_CHECK( phrase_parse(beg, end, defp, space, result) );
  BOOST_CHECK( (beg == end) );                         // we should consume all input
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
  std::string testdef("DESIGN test ;\nCOMPONENTS 1 ;\n - I111 INVX2 ;\nEND COMPONENTS\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  def result;
  BOOST_CHECK( phrase_parse(beg, end, defp, space, result) );
  BOOST_CHECK( (beg == end) );                         // we should consume all input
  BOOST_CHECK_EQUAL( result.name, "test" );
  BOOST_REQUIRE_EQUAL( result.components.size(), 1 );    // exactly one component read
  BOOST_CHECK_EQUAL( result.components[0].name, "I111" );
  BOOST_CHECK_EQUAL( result.components[0].celltype, "INVX2" );
  BOOST_CHECK( (!result.components[0].placement) );     // it has NO placement
}  

BOOST_AUTO_TEST_CASE ( components_parse_wrongcount ) {
  std::string testdef("DESIGN test ;\nCOMPONENTS 2 ;\n - I111 INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  BOOST_CHECK( !phrase_parse(beg, end, defp, space) );
}
  
// BOZO when we eventually parse everything this won't be a very interesting test and probably should be removed,
// or have previously ignored stuff checked
BOOST_AUTO_TEST_CASE ( parse_ignored_stuff ) {
  std::string testdef("DESIGN test ;\nVERSION 1.211 ;\nDIEAREA ( 0 0 ) ( 100000 200000 ) ;\nCOMPONENTS 1 ;\n - I111 INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nSITE CORE1 0 0 N DO 200 BY 1 STEP 100 500 ;\nSPECIALNETS 1 ;\n - GND ;\nEND SPECIALNETS\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  def result;
  BOOST_CHECK( phrase_parse(beg, end, defp, space, result) );
  BOOST_CHECK( (beg == end) );                         // we should consume all input
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
  
