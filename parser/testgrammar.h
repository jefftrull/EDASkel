#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

// a starter grammar
template <typename Iterator>
struct testparser : boost::spirit::qi::grammar<Iterator,
                                              boost::spirit::qi::space_type>
{

 testparser() : testparser::base_type(elt_section)
    {
      using namespace boost::spirit::qi;
      using boost::spirit::_1;                  // access attributes for component count check
      using boost::spirit::_r1;
      using boost::spirit::_r2;
      using boost::phoenix::val;                // for error handling
      using boost::phoenix::construct;          // for error handling
      using boost::spirit::qi::locals;
      int n;                                    // for counting
      using boost::phoenix::ref;                // for counting

      // the repeated items
      element = '-' > lexeme[alpha > *alnum] > ';' ;

      // a specific rule
      specific_count_elements = lit("ELTTYPE") > int_[_a = _1] > ';' >  // remember count
	repeat(_a)[element] >                        // expect that many
	lit("END") > lit("ELTTYPE") ;

      // a general rule
      count_elements = lit(_r1) > int_[_a = _1] > ';' >  // remember count
                  	repeat(_a)[_r2] >                        // expect that many
                	lit("END") > lit(_r1) ;

      // instances of the general rule
      // COMPILE ERROR HERE
      elt_section = count_elements(val("ELTTYPE"), ref(element)) ;


    }

  typedef boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<int>, boost::spirit::qi::space_type> LocalRule;
  typedef boost::spirit::qi::rule<Iterator, boost::spirit::qi::space_type> Rule;
  typedef boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<int>, boost::spirit::qi::space_type, void(std::string, Rule) > CountRule;
  LocalRule specific_count_elements;
  Rule element, elt_section;
  CountRule count_elements;
};

