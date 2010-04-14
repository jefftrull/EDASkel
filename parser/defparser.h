// A Spirit-based DEF parser
// Copyright 2010 Jeffrey Elliot Trull <linmodemstudent@gmail.com>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include "keyword.h"
#include "deftypes.h"

#include <vector>

// a starter DEF grammar
template <typename Iterator>
struct defparser : boost::spirit::qi::grammar<Iterator,
                                              std::vector<defcomponent>(),
                                              boost::spirit::qi::locals<int>, 
                                              boost::spirit::qi::space_type>
{

  defparser() : defparser::base_type(def_file)
    {
      using namespace boost::spirit::qi;
      using namespace distinct;
      using boost::spirit::_1;                  // access attributes for component count check
      using boost::phoenix::val;                // for error handling
      using boost::phoenix::construct;          // for error handling
      using boost::phoenix::ref;                // for counting

      version_stmt = keyword["VERSION"] >> double_ >> ';' ;
      diearea_stmt = keyword["DIEAREA"] >> '(' >> int_ >> int_ >> ')' >> '(' >> int_ >> int_ >> ')' ;

      // define some major components

      // components (instances)
      orient %= lexeme[- char_("F") >> char_("NSEW")] ;
      // Using ">" here instead of ">>" implies a required sequence and allows the parser to check it
      // specifically (instead of simply failing on the whole input)
      plcinfo %= '+' > (string("FIXED") | string("PLACED")) >
		     '(' > int_ > int_ > ')' > orient > ';' ;    // location and orientation

      component %= '-' > lexeme[alpha > *alnum] > lexeme[alpha > *alnum] >
	           - plcinfo ;

      // TBD rule - absorb everything up to the trailing semicolon
      tbd_elt = '-' > *char_ > ';' ;

      // define repeat rules for each element section
      // I tried to make this generic but failed.  You can pass in the rule as an inherited attribute,
      // but the type of the parent rule needs to be different for each since the synthesized attribute is
      // different (i.e., std::vector<attribute_of_child_rule>)
      comps_section %= keyword["COMPONENTS"] > omit[int_[_a = _1]] > ';' >  // remember count in a local variable
	               repeat(_a)[component] >                        // expect that many copies of the supplied rule
	               keyword["END"] > keyword["COMPONENTS"] ;       // END followed by the section name again

      pins_section = keyword["PINS"] > int_[_a = _1] > ';' >
	             repeat(_a)[tbd_elt] >
	             keyword["END"] > keyword["PINS"] ;

      vias_section = keyword["VIAS"] > int_[_a = _1] > ';' >
	             repeat(_a)[tbd_elt] >
	             keyword["END"] > keyword["VIAS"] ;

      nets_section = keyword["NETS"] > int_[_a = _1] > ';' >
	             repeat(_a)[tbd_elt] >
	             keyword["END"] > keyword["NETS"] ;

      snets_section = keyword["SPECIALNETS"] > int_[_a = _1] > ';' >
	              repeat(_a)[tbd_elt] >
	              keyword["END"] > keyword["SPECIALNETS"] ;

      // a generic rule describing how the multi-element sections (COMPONENTS, PINS, VIAS, etc.) work
      // this is a rule with two "inherited parameters" (i.e., arguments), one of which is the name
      // of the section, and one of which is the rule describing the element
      // this doesn't work, actually, because we have to define the synthesized attribute differently for each caller
      counted_elements %= keyword[_r1] > int_[_a = _1] > ';' >  // remember count in a local variable
	                 repeat(_a)[_r2] >                     // expect that many copies of the supplied rule
	                 keyword["END"] > keyword[_r1] ;       // END followed by the section name again

      // apply the generic rule to several examples
      //      comps_section %= counted_elements(val("COMPONENTS"), ref(component)) ;

      // of the top-level statements, some can occur only once (VERSION) and others can occur many times (SITE)
      // if not for that we could use permutation operator (^) which would check uniqueness for us
      // topstmt %= version_stmt | comps_section | nets_section | diearea_stmt ;
      topstmt %= version_stmt | comps_section | nets_section | diearea_stmt ;

      def_file %= keyword["DESIGN"] >> omit[lexeme[+alnum]] >> ';' >> topstmt >> keyword["END"] >> keyword["DESIGN"] ;
      // def_file %= keyword["DESIGN"] >> omit[lexeme[+alnum]] >> ';' >> *topstmt >> keyword["END"] >> keyword["DESIGN"] ;
      // BOZO check there is at most one of components, nets, specialnets, pins, version, etc.

        on_error<fail>
        (
	 def_file
	 , std::cerr
	 << val("Error! Expecting ")
	 << boost::spirit::_4                               // what failed?
	 << val(" here: \"")
	 << construct<std::string>(boost::spirit::_3, boost::spirit::_2)   // iterators to error-pos, end
	 << val("\"")
	 << std::endl
        );
    }

  typedef boost::spirit::qi::rule<Iterator, boost::spirit::qi::space_type> Rule;
  Rule end_design;
  Rule version_stmt, diearea_stmt;
  Rule tbd_elt;
  typedef boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::qi::space_type> StringRule;
  StringRule orient;
  typedef boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<int>, boost::spirit::qi::space_type> LocalVarRule;
  LocalVarRule pins_section, nets_section, snets_section, vias_section;

  // optional placement info for a component
  typedef boost::spirit::qi::rule<Iterator, defplcinfo(), boost::spirit::qi::space_type > PlcRule;
  PlcRule plcinfo;

  // a component rule produces a defcomponent struct as its synthesized attribute
  typedef boost::spirit::qi::rule<Iterator, defcomponent(), boost::spirit::qi::space_type > CompRule;
  CompRule component;

  // a component rule producing a vector of defcomponents
  typedef boost::spirit::qi::rule<Iterator, std::vector<defcomponent>(), boost::spirit::qi::locals<int>, boost::spirit::qi::space_type > CompVecRule;
  CompVecRule comps_section, topstmt, def_file;   // BOZO def_file is more general, but leave it for now

  typedef boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<int>, void(std::string, Rule), boost::spirit::qi::space_type > CountRule;
  CountRule counted_elements;
};

