#define BOOST_TEST_MODULE basic test
#include <boost/test/included/unit_test.hpp>

#include "../defparser.h"
#include <string>
#include <iostream>

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
  BOOST_CHECK( phrase_parse(beg, end, defp, space) );
  BOOST_CHECK( (beg == end) );                         // we should consume all input

  // BOZO verify 0 components stored
}
  
BOOST_AUTO_TEST_CASE ( components_parse_simple ) {
  std::string testdef("DESIGN test ;\nCOMPONENTS 1 ;\n - I111 INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  std::vector<defcomponent> result;
  BOOST_CHECK( phrase_parse(beg, end, defp, space, result) );
  BOOST_CHECK( (beg == end) );                         // we should consume all input
  BOOST_CHECK_EQUAL( result.size(), 1 );
  // BOZO verify location, orientation, placement type, name of component
}
  
BOOST_AUTO_TEST_CASE ( components_parse_wrongcount ) {
  std::string testdef("DESIGN test ;\nCOMPONENTS 2 ;\n - I111 INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  BOOST_CHECK( !phrase_parse(beg, end, defp, space) );
}
  
// BOZO when we parse everything this test goes away
BOOST_AUTO_TEST_CASE ( parse_ignored_stuff ) {
  std::string testdef("DESIGN test ;\nVERSION 1.211 ;\nDIEAREA ( 0 0 ) ( 100000 200000 ) ;\nCOMPONENTS 1 ;\n - I111 INVX2 + FIXED ( -4107 82000 ) FN ;\nEND COMPONENTS\nSITE CORE1 0 0 N DO 200 BY 1 STEP 100 500 ;\nEND DESIGN\n");
  std::string::const_iterator beg = testdef.begin();
  std::string::const_iterator end = testdef.end();
  BOOST_CHECK( phrase_parse(beg, end, defp, space) );
  BOOST_CHECK( (beg == end) );                         // we should consume all input
  // BOZO check that version, diearea, and one component are present
}
  
