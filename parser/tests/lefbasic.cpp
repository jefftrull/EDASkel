// Tests for a Boost Spirit-based LEF parser, part of EDASkel, a sample EDA app
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
#define BOOST_TEST_MODULE basic LEF tests
#include <boost/test/included/unit_test.hpp>

#include "../lefparser.h"
#include <string>
#include <iostream>

lefparser<std::string::const_iterator> lefp;
lefskipper<std::string::const_iterator> skp;

BOOST_AUTO_TEST_CASE( case_check ) {

  std::string testlef("NAMESCASESENSITIVE ON ;");
  std::string::const_iterator beg = testlef.begin();
  std::string::const_iterator end = testlef.end();
  lef result;
  BOOST_CHECK( phrase_parse(beg, end, lefp, skp, result) );  // we should match
  BOOST_CHECK( (beg == end) );                         // we should consume all input
  // did not use BOOST_CHECK_EQUAL b/c it wants to output these on failure, and there is no operator<< defined
  BOOST_CHECK( result.macros.empty() );

}

BOOST_AUTO_TEST_CASE( macro_basic_check ) {

  std::string testlef("MACRO INX2\nCLASS CORE ;\n FOREIGN INX2 0.0 -1.0 ;\nORIGIN 0.0 1.0 ;\nSIZE 1.0 BY 10.0 ;\nSYMMETRY X Y ;\nEND INX2");
  std::string::const_iterator beg = testlef.begin();
  std::string::const_iterator end = testlef.end();
  lef result;
  BOOST_CHECK( phrase_parse(beg, end, lefp, skp, result) );  // we should match
  BOOST_CHECK( (beg == end) );                         // we should consume all input
  if (beg != end)
    std::cerr << "remaining input data is as follows:|" << std::string(beg, end) << "|\n";
  // did not use BOOST_CHECK_EQUAL b/c it wants to output these on failure, and there is no operator<< defined
  BOOST_REQUIRE_EQUAL( result.macros.size(), 1 );
  BOOST_CHECK_EQUAL( result.macros[0].name, "INX2" );
  BOOST_CHECK_EQUAL( result.macros[0].class_, SITECLASS_CORE );
  BOOST_REQUIRE( result.macros[0].foreign );
  BOOST_CHECK_EQUAL( result.macros[0].foreign->name, "INX2" );
  BOOST_CHECK_CLOSE( result.macros[0].foreign->pt.x, 0.0, 0.001f );
  BOOST_CHECK_CLOSE( result.macros[0].foreign->pt.y, -1.0, 0.001f );
  BOOST_REQUIRE( result.macros[0].origin );
  BOOST_CHECK_CLOSE( result.macros[0].origin->x, 0.0, 0.001f );
  BOOST_CHECK_CLOSE( result.macros[0].origin->y, 1.0, 0.001f );
  BOOST_REQUIRE( result.macros[0].size );
  BOOST_CHECK_CLOSE( result.macros[0].size->width, 1.0, 0.001f );
  BOOST_CHECK_CLOSE( result.macros[0].size->height, 10.0, 0.001f );
  BOOST_REQUIRE( result.macros[0].symmetry );
  BOOST_REQUIRE_EQUAL( result.macros[0].symmetry->size(), 2 );
  // BOOST_CHECK_EQUAL( result.macros[0].symmetry->[0], SITESYM_X );
  // BOOST_CHECK_EQUAL( result.macros[0].symmetry->[1], SITESYM_Y );
}

BOOST_AUTO_TEST_CASE( site_basic_check ) {

  std::string testlef("SITE MYSITENAME CLASS PAD ; SYMMETRY R90 ; SIZE 11.01 BY 22 ; END MYSITENAME");
  std::string::const_iterator beg = testlef.begin();
  std::string::const_iterator end = testlef.end();
  lef result;
  BOOST_CHECK( phrase_parse(beg, end, lefp, skp, result) );
  BOOST_CHECK( (beg == end) );
  BOOST_REQUIRE_EQUAL( result.sites.size(), 1 );
  BOOST_CHECK_EQUAL( result.sites[0].name, "MYSITENAME" );
  BOOST_CHECK_EQUAL( result.sites[0].class_, SITECLASS_PAD );
  BOOST_REQUIRE( result.sites[0].symmetry );
  BOOST_CHECK_EQUAL( result.sites[0].symmetry, SITESYM_R90 );
  BOOST_CHECK_CLOSE( result.sites[0].width, 11.01, 0.001f );
  BOOST_CHECK_CLOSE( result.sites[0].height, 22.0, 0.001f );

}

BOOST_AUTO_TEST_CASE( reproduce_quasiHP ) {

  std::string testlef("#\n# CMOS26G LEF technology file\nNAMESCASESENSITIVE ON ;\n#UNITS\n  # not allowed in lef: DISTANCE MICRON 1 ;\nLAYER poly\n# no routing will be done in poly\n  TYPE MASTERSLICE ;\nEND poly\n\nVIA via01 DEFAULT\n  #eas FOREIGN VIA01 0 0 ;\n  LAYER poly      ; RECT -0.8 -0.8 0.8 0.8 ; # 8A1 & 8G\n  LAYER cont1     ; RECT -0.4 -0.4 0.4 0.4 ; # 8A1\n  LAYER metal1    ; RECT -0.8 -0.8 0.8 0.8 ; # 8A1 & 8E1\n  RESISTANCE 8.0 ; # nom at 110C\nEND via01\n\nVIARULE pwrstripe\n  layer metal1 ; direction horizontal ; width 5.6 to 5.6 ;\n  layer metal2 ; direction vertical ; width 6.8 to 6.8 ;\n  via stripevia ;\nEND pwrstripe\n\nSPACING\n  SAMENET cont2 cont2 1.2 ; # Rule 10B\n  SAMENET cont1 cont2 0.6 ; # Rule 10C1\n  SAMENET cont3 cont3 1.2 ; # Rule 12B\n  SAMENET cont2 cont3 0.6 ; # Rule 12C\nEND SPACING\n\nSITE CORE26  CLASS CORE ; SYMMETRY    Y ; SIZE   2.6 BY  38.4 ; END CORE26\n#Site for dummy IO ports\nSITE SMALLIO  CLASS PAD ;                SIZE   1.3 BY   1.3 ; END SMALLIO\n\nMACRO AND2EE\n  CLASS CORE ;\n  FOREIGN AND2EE 0.000 -4.800 ;\n  ORIGIN 0.000 4.800 ;\n  SIZE 10.400 BY 38.400 ;\n  SYMMETRY X Y ;\n  SITE CORE26 ;\n  PIN A\n    DIRECTION INPUT ;\n    USE SIGNAL ;\n    CAPACITANCE 0.00 ;\n    PORT\n    LAYER metal1 ;\n      RECT 8.600 9.100 9.600 12.500 ; # A|0.0@0\n    END\n  END A\n  PIN B\n    DIRECTION INPUT ;\n    USE SIGNAL ;\n    CAPACITANCE 0.00 ;\n    PORT\n    LAYER metal1 ;\n      RECT 6.000 9.100 7.000 12.500 ; # B|0.0@0\n    END\n  END B\n  PIN Q\n    DIRECTION OUTPUT ;\n    USE SIGNAL ;\n    CAPACITANCE 0.00 ;\n    PORT\n    LAYER metal1 ;\n      RECT 0.800 6.700 1.800 19.700 ; # Q|0.0@0\n    END\n  END Q\n  PIN VDD\n    DIRECTION INOUT ;\n    USE POWER ;\n    SHAPE ABUTMENT ;\n    CAPACITANCE 0.00 ;\n    PORT\n    LAYER metal1 ;\n      RECT 0.000 23.600 10.400 29.200 ; # VDD|0.0@0\n    END\n  END VDD\n  PIN GND\n    DIRECTION INOUT ;\n    USE GROUND ;\n    SHAPE ABUTMENT ;\n    CAPACITANCE 0.00 ;\n    PORT\n    LAYER metal1 ;\n      RECT 0.000 -0.400 10.400 5.200 ; # GND|0.0@0\n    END\n  END GND\n  OBS\n    LAYER metal1 ;\n      RECT 3.400 6.700 9.600 7.700 ;\n      RECT 3.400 6.700 4.400 17.300 ;\n      RECT 3.400 16.300 7.000 17.300 ;\n      RECT 6.000 16.300 7.000 19.700 ;\n    LAYER cont2 ;\n      RECT 3.400 -0.500 7.000 0.500 ;\n      RECT 0.800 4.300 1.800 7.700 ;\n      RECT 8.600 4.300 9.600 7.700 ;\n      RECT 0.800 6.700 9.600 7.700 ;\n      RECT 3.400 6.700 4.400 19.700 ;\n      RECT 3.400 11.500 9.600 12.500 ;\n      RECT 3.400 16.300 7.000 19.700 ;\n      RECT 0.800 18.700 7.000 19.700 ;\n      RECT 3.400 23.500 9.600 24.500 ;\n      RECT 3.400 28.300 4.400 29.300 ;\n  END\nEND AND2EE\n\nEND LIBRARY\n# /mnts/hpdtctm/cbs1/D3.1X/700/bin/gridit:\n# 	 set-b.ada 1.6 93/05/05 15:24:33\n");
  std::string::const_iterator beg = testlef.begin();
  std::string::const_iterator end = testlef.end();
  lef result;
  BOOST_CHECK( phrase_parse(beg, end, lefp, skp, result) );
  BOOST_CHECK( (beg == end) );
  if (beg != end)
    std::cerr << "remaining input data is as follows:|" << std::string(beg, end) << "|\n";
}
