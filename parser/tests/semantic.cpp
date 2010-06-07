// Tests for semantic analysis for Boost Spirit-based LEF/DEF parser, part of EDASkel, a sample EDA app
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

#define BOOST_TEST_MODULE basic test
#include <boost/test/included/unit_test.hpp>

#include "../../db/simpledb.h"
#include "../defparser.h"
#include "../defsem.h"
#include <string>
#include <iostream>

using namespace DefParse;
using namespace SimpleDB;

defparser<std::string::const_iterator> defp;
using boost::spirit::qi::space;

BOOST_AUTO_TEST_CASE ( diearea_checks ) {
  // parse syntax into "def" parse syntax structure
  std::string testdef("DESIGN test ;\nDIEAREA ( -2000 -2000 ) ( 100000 200000 ) ;\nEND DESIGN");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  def result;
  BOOST_CHECK( phrase_parse(beg, end, defp, space, result) );

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
  testdef = std::string("DESIGN test ;\nDIEAREA ( 10 10 ) ( 0 0 ) ;\nEND DESIGN");
  beg = testdef.begin(); end = testdef.end();
  def result_badboundary;
  BOOST_CHECK( phrase_parse(beg, end, defp, space, result_badboundary) );

  Database db_badboundary;
  // our default policy says we don't abort but we do skip errors
  BOOST_CHECK( chk.CheckAndInsert(result_badboundary, lib, db_badboundary) ); // no abort
  BOOST_CHECK( !db_badboundary.hasExtent() );              // boundary NOT set

}

// change the policy (by specializing for diearea) and verify it takes effect
template<DefCheckError Err>
class DieAreaAbortPolicy : public DefCheckPolicy<Err> {};

template<>
class DieAreaAbortPolicy<DEFERR_DIE_MALFORMED_RECT> {
public:
    static bool const silent = true;  // one less message during test runs
    static bool const skip = false;   // actually set the rect, BUT
    static bool const abort = true;   // also DO abort
};

BOOST_AUTO_TEST_CASE ( policy_checks ) {
  std::string testdef("DESIGN test ;\nDIEAREA ( 10 10 ) ( 0 0 ) ;\nEND DESIGN");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  def result;
  BOOST_CHECK( phrase_parse(beg, end, defp, space, result) );

  // turn syntax result into database contents (while checking)
  Library lib;
  Database db;
  DefChecker<Database, Library, DieAreaAbortPolicy> chk;
  BOOST_CHECK( !chk.CheckAndInsert(result, lib, db) ); // check should abort
  BOOST_CHECK( db.hasExtent() );                   // but boundary should be present
}
