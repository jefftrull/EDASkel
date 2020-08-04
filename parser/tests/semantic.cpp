// Tests for semantic analysis for Boost Spirit-based LEF/DEF parser, part of EDASkel, a sample EDA app
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

#define BOOST_TEST_MODULE basic test
#include <boost/test/included/unit_test.hpp>

#include <string>
#include <iostream>

#include "../../db/simpledb.hpp"
#include "../lefparser.hpp"
#include "../lefsem.hpp"
#include "../defparser.hpp"
#include "../defsem.hpp"
#include "../lefdef.h"

using namespace LefParse;
using namespace DefParse;
using namespace SimpleDB;

namespace EDASkel {
  extern DefTokens<LefDefLexer> defTokens;
  extern defparser<DefTokens<LefDefLexer>::iterator_type, DefTokens<LefDefLexer>::lexer_def > defParser;
  extern LefTokens<LefDefLexer> lefTokens;
  extern lefparser<LefTokens<LefDefLexer>::iterator_type, LefTokens<LefDefLexer>::lexer_def > lefParser;
}
using namespace boost::spirit::qi;

BOOST_AUTO_TEST_CASE ( diearea_checks ) {
  // parse syntax into "def" parse syntax structure
  std::stringstream testdef("DESIGN test ;\nDIEAREA ( -2000 -2000 ) ( 100000 200000 ) ;\nEND DESIGN");
  testdef.unsetf(std::ios::skipws);
  LefDefIter beg = LefDefIter(testdef), end;
  def result;
  BOOST_CHECK( tokenize_and_parse(beg, end, defTokens, defParser, result) );

  // turn syntax result into database contents (while checking)
  Library lib;
  Database db;
  DefChecker<Database, Library> chk;
  BOOST_CHECK( chk.CheckAndInsert(result, lib, db) );

  // verify database contents
  Database::Rect boundary = db.getExtent();
  BOOST_CHECK( (boundary.ll().x() == -2000) && (boundary.ll().y() == -2000) &&
	       (boundary.ur().x() == 100000) && (boundary.ur().y() == 200000) );

  // BOZO add instance and site checks

  // now try a bad diearea and verify the check fails
  std::stringstream testdefbad("DESIGN test ;\nDIEAREA ( 10 10 ) ( 0 0 ) ;\nEND DESIGN");
  testdefbad.unsetf(std::ios::skipws);
  beg = LefDefIter(testdefbad);
  def result_badboundary;
  BOOST_CHECK( tokenize_and_parse(beg, end, defTokens, defParser, result_badboundary) );

  Database db_badboundary;
  // our default policy says we don't abort but we do skip errors
  BOOST_CHECK( chk.CheckAndInsert(result_badboundary, lib, db_badboundary) ); // no abort
  BOOST_CHECK( !db_badboundary.hasExtent() );              // boundary NOT set

}

// change the policy (by specializing for diearea) and verify it takes effect
template<DefCheckError Err>
class DieAreaAbortPolicy : public DefCheckPolicy<Err> {};

template<>
class DieAreaAbortPolicy<DEFERR_DIEAREA_MALFORMED> {
public:
    static bool const silent = true;  // one less message during test runs
    static bool const skip = false;   // actually set the rect, BUT
    static bool const abort = true;   // also DO abort
};

BOOST_AUTO_TEST_CASE ( policy_checks ) {
  std::stringstream testdef("DESIGN test ;\nDIEAREA ( 10 10 ) ( 0 0 ) ;\nEND DESIGN");
  testdef.unsetf(std::ios::skipws);
  LefDefIter beg = LefDefIter(testdef), end;
  def result;
  BOOST_CHECK( tokenize_and_parse(beg, end, defTokens, defParser, result) );

  // turn syntax result into database contents (while checking)
  Library lib;
  Database db;
  DefChecker<Database, Library, DieAreaAbortPolicy> chk;
  BOOST_CHECK( !chk.CheckAndInsert(result, lib, db) ); // check should abort
  BOOST_CHECK( db.hasExtent() );                   // but boundary should be present
}

BOOST_AUTO_TEST_CASE ( lefdef_combined_basic ) {
  // LEF: one site and one macro
  std::stringstream testlef("SITE CORE0 CLASS CORE ; SYMMETRY Y ; SIZE 1.0 BY 2.0 ; END CORE0\nMACRO INX2\nCLASS CORE ;\nORIGIN 0.0 0.0 ;\nSIZE 2.0 BY 2.0 ;\nSYMMETRY X Y ;\nSITE CORE0 ;\nEND INX2");
  testlef.unsetf(std::ios::skipws);
  LefDefIter beg(testlef), end;
  lef lefresult;
  BOOST_CHECK( tokenize_and_parse(beg, end, lefTokens, lefParser, lefresult) );
  BOOST_CHECK( beg == end );
  Library lib;
  LefChecker<Library> lchk;
  BOOST_CHECK( lchk.CheckAndInsert(lefresult, lib) );
  // specific library checks here
  using CellPtr = Library::CellPtr;
  CellPtr cell = lib.findCell("INX2");
  BOOST_REQUIRE( cell != CellPtr() );
  BOOST_CHECK_CLOSE( cell->getWidth(), 2.0, 0.001f );
  BOOST_CHECK_CLOSE( cell->getHeight(), 2.0, 0.001f );
  BOOST_CHECK( cell->getSymmetry().size() == 2 );

  using SitePtr = Library::SitePtr;
  SitePtr site = lib.findSite("CORE0");
  BOOST_REQUIRE( site != SitePtr() );
  BOOST_REQUIRE( site->getSymmetry().size() == 1 );
  BOOST_CHECK( site->getSymmetry()[0] == SITESYM_Y );
  BOOST_CHECK_CLOSE( site->getWidth(), 1.0, 0.001f );
  BOOST_CHECK_CLOSE( site->getHeight(), 2.0, 0.001f );

  // DEF: define a diearea with set of sites and put two macros there
  std::stringstream testdef("DESIGN test ;\nDIEAREA ( 0 0 ) ( 1000 1000 ) ;\nUNITS DISTANCE MICRONS 100 ;\nSITE CORE0 0 0 N DO 10 BY 1 STEP 100 200 ;\nSITE CORE0 0 200 FS DO 10 BY 1 STEP 100 200 ;\nSITE CORE0 0 400 N DO 10 BY 1 STEP 100 200 ;\nSITE CORE0 0 600 FS DO 10 BY 1 STEP 100 200 ;\nSITE CORE0 0 800 N DO 10 BY 1 STEP 100 200 ;\nCOMPONENTS 2 ;\n- inst1 INX2 + PLACED ( 300 600 ) FS ;\n- inst2 INX2 + PLACED ( 500 800 ) N ;\nEND COMPONENTS\nEND DESIGN");
  testdef.unsetf(std::ios::skipws);
  def defresult;
  beg = LefDefIter(testdef);
  BOOST_CHECK( tokenize_and_parse(beg, end, defTokens, defParser, defresult) );
  BOOST_CHECK( beg == end );

  DefChecker<Database, Library, DieAreaAbortPolicy> chk;
  Database db;
  BOOST_CHECK( chk.CheckAndInsert(defresult, lib, db) );

  // verify location, orientation, etc. of sites and placed instances
  boost::optional<std::string> sitename = db.siteAt(0, 0);
  BOOST_REQUIRE( sitename );
  BOOST_CHECK( *sitename == "CORE0" );
  sitename = db.siteAt(100, 200);
  BOOST_REQUIRE( sitename );
  BOOST_CHECK( *sitename == "CORE0" );
  // extreme upper right site
  sitename = db.siteAt(900, 800);
  BOOST_REQUIRE( sitename );
  BOOST_CHECK( *sitename == "CORE0" );
  // just outside - bottom right
  sitename = db.siteAt(1000, 0);
  BOOST_CHECK( !sitename );
  // just outside - upper left
  sitename = db.siteAt(0, 1000);
  BOOST_CHECK( !sitename );
  // TBD: is it valid to have a cell not on a site boundary, if it completely overlaps sites?
  // Our existing logic assumes this is *not* valid
  // off-grid in X:
  sitename = db.siteAt(50, 200);
  BOOST_CHECK( !sitename );
  // off-grid in Y:
  sitename = db.siteAt(100, 300);
  BOOST_CHECK( !sitename );

  using InstPtr = Database::InstPtr;
  InstPtr inst = db.findInst("inst1");
  BOOST_REQUIRE( inst != InstPtr() );
  BOOST_CHECK( inst->getCellName() == "INX2" );
  BOOST_REQUIRE( inst->hasPlacement() );
  BOOST_CHECK( (inst->getOrigin().x() == 300) && (inst->getOrigin().y() == 600) );
  BOOST_CHECK( inst->getOrient() == "FS" );
  BOOST_CHECK( !inst->isFixed() );
  sitename = db.siteAt(300, 600);
  BOOST_REQUIRE( sitename );
  BOOST_CHECK( *sitename == "CORE0" );

  inst = db.findInst("inst2");
  BOOST_REQUIRE( inst != InstPtr() );
  BOOST_CHECK( inst->getCellName() == "INX2" );
  BOOST_REQUIRE( inst->hasPlacement() );
  BOOST_CHECK( (inst->getOrigin().x() == 500) && (inst->getOrigin().y() == 800) );
  BOOST_CHECK( inst->getOrient() == "N" );
  BOOST_CHECK( !inst->isFixed() );
}
