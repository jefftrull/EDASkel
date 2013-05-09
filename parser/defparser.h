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
#if !defined(EDASKEL_DEF_PARSER)
#define EDASKEL_DEF_PARSER

#include <vector>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <boost/spirit/repository/include/qi_kwd.hpp>
#include <boost/spirit/repository/include/qi_keywords.hpp>

#include "deftypes.h"
#include "lefdef.h"

using namespace EDASkel;

namespace DefParse {

// a starter DEF grammar

// first define some major subcomponents in their own "grammar"

// Components
typedef boost::spirit::qi::symbols<char, std::string> comp_symtab_t;  // for reference checking

template <typename Iterator>
struct comp_parser : boost::spirit::qi::grammar<Iterator,
                                              defcomponent(),
                                              lefdefskipper<Iterator> >
{

  comp_parser(comp_symtab_t& compsym) : comp_parser::base_type(component), comp_symtab(compsym)
    {
      using namespace boost::spirit::qi;
      using boost::spirit::repository::dkwd;
      using boost::phoenix::val;                // for error handling
      using boost::phoenix::construct;          // for error handling
      using boost::phoenix::at_c;
      using boost::phoenix::bind;
      using boost::spirit::_1;
      using boost::spirit::_2;
      using boost::spirit::_val;

      point %= '(' >> int_ >> int_ >> ')' ;     // points are parenthesized pairs, no comma

      // components (instances)
      orient = repeat(1)[char_("NSEW")] |       // a single cardinal direction, OR
	       (char_('F') >> char_("NSEW")) ;  // the same, but flipped

      // Using ">" here instead of ">>" implies a required sequence and allows the parser to check it
      // specifically (instead of simply failing on the whole input)
      plcinfo %= '+' >> (dkwd("FIXED", 1)[attr("FIXED") > point > orient] |
			 dkwd("PLACED", 1)[attr("PLACED") > point > orient]) ;    // location and orientation

      weight %= '+' >> dkwd("WEIGHT", 1)[int_] ;

      source = '+' >> dkwd("SOURCE", 1) [(dkwd("DIST")[attr("DIST")] /
					  dkwd("NETLIST")[attr("NETLIST")] /
					  dkwd("USER")[attr("USER")] /
					  dkwd("TIMING")[attr("TIMING")])] ;

      // instance names: letter followed by letters, numbers, underscore, hyphen, square brackets, slashes (hierarchy)
      // BOZO support escapes
      iname %= lexeme[alpha >> *(alnum | char_('-') | char_('_') | char_('[') | char_(']') | char_('/'))] ;
      // celltypes: assuming only underscore might be used out of the non-alphanumeric chars
      ctype %= lexeme[alpha >> *(alnum | char_('_'))] ;

      // components required instance name and celltype; optional any (or none) of placement and weight, in any order:
      component = ('-' > iname[at_c<0>(_val) = _1] > ctype[at_c<1>(_val) = _1] >
                   (plcinfo[at_c<2>(_val) = _1] ^ omit[weight] ^ source ^ eps ) >
                   ';' )[bind(comp_symtab.add, _1, _1)];  // turn symbol table add into "lazy" function

      weight.name("Weight");
      source.name("Source");
      iname.name("Instance Name");
      ctype.name("Cell Type");
      plcinfo.name("Optional Placement Info");
      component.name("Component");

      on_error<fail>
        (
	 component
       , std::cerr
       << val("Error! Expecting ")
       << boost::spirit::_4                               // what failed?
       << val(" here: \"")
       << construct<std::string>(boost::spirit::_3, boost::spirit::_2)   // iterators to error-pos, end
       << val("\"")
       << std::endl
        );

      // BOZO testing
      comp_symtab.add("foo", "foo");
      comp_symtab.add(std::string("bar"), std::string("bar"));

    }

  // helpful abbreviations
  typedef lefdefskipper<Iterator> skipper;

  template<typename Signature>
    struct Rule
    {
      typedef boost::spirit::qi::rule<Iterator, Signature, skipper> type;
    };

  // boost::spirit::qi::rules pertaining to COMPONENTS - see deftypes.h for data results

  // points "( x y )" produces defpoint structs (see deftypes.h)
  typename Rule<defpoint()>::type point;

  // optional placement info (placed vs. fixed, location, orientation)
  typename Rule<defplcinfo()>::type plcinfo;

  typename Rule<std::string()>::type orient, iname, ctype;


  // boost::spirit::qi::rules neither inheriting nor synthesizing an attribute
  // BOZO try unused_type here
  boost::spirit::qi::rule<Iterator, skipper> weight, source;

  // a single instance within the COMPONENTS section (name, celltype, placement)
  typename Rule<defcomponent()>::type component;

  // Symbol table storage for components
  comp_symtab_t& comp_symtab;

};

// A NET statement
template <typename Iterator>
struct net_parser : boost::spirit::qi::grammar<Iterator,
                                              defnet(),
                                              lefdefskipper<Iterator> >
{
  net_parser(comp_symtab_t const& compsym) : net_parser::base_type(net), comp_symtab(compsym)
   {
     using namespace boost::spirit::qi;

     nname %= lexeme[alpha >> *(alnum | char_('-') | char_('_') | char_('[') | char_(']') | char_('/'))] ;

     connection = '(' > comp_symtab > nname > ')' ;

     net = '-' > nname > *connection > ';' ;
  }

  typedef lefdefskipper<Iterator> skipper_t;
  boost::spirit::qi::rule<Iterator, std::string(), skipper_t> nname;
  boost::spirit::qi::rule<Iterator, std::pair<std::string, std::string>(), skipper_t> connection;
  boost::spirit::qi::rule<Iterator, defnet(), skipper_t> net;
  comp_symtab_t const& comp_symtab;
};  


template <typename Iterator>
struct defparser : boost::spirit::qi::grammar<Iterator,
                                              boost::spirit::qi::locals<int, std::string>,
                                              def(),
                                              lefdefskipper<Iterator> >
{

  defparser() : defparser::base_type(def_file), component(comp_symtab), net(comp_symtab)
    {
      using namespace boost::spirit::qi;

      using boost::spirit::_1;                  // access attributes for component count check
      using boost::phoenix::push_back;          // to store results in containers
      using boost::phoenix::val;                // for error handling
      using boost::phoenix::construct;          // for error handling
      using boost::phoenix::at_c;               // to refer to pieces of wrapped structs
      using boost::spirit::repository::dkwd;
      using boost::spirit::repository::kwd;
      using boost::spirit::qi::no_skip;         // when we (rarely) actually want the whitespace

      // top-level elements in a DEF file
      version_stmt = dkwd("VERSION", 1)[double_] > ';' ;

      point %= '(' >> int_ >> int_ >> ')' ;     // points are parenthesized pairs, no comma
      rect %= point >> point ;                  // rects are just two points in a row
      diearea_stmt %= dkwd("DIEAREA", 1)[rect] > ';' ;

      // identifiers
      // design name: letter followed by letters, numbers, underscore, or hyphen
      dname %= lexeme[alpha >> *(alnum | char_('-') | char_('_'))] ;

      // celltypes: assuming only underscore might be used out of the non-alphanumeric chars
      ctype %= lexeme[alpha >> *(alnum | char_('_'))] ;

      // orientation
      orient = repeat(1)[char_("NSEW")] |       // a single cardinal direction, OR
	        (char_('F') >> char_("NSEW")) ; // the same, but flipped

      // define some major elements

      // My copy of the LEF/DEF reference does not show this SITE command as valid for DEF yet my example data does...
      // The example data's syntax is very similar to that defined for ROW, so I'll combine them
      siterpt_stmt = dkwd("DO")[int_] > dkwd("BY")[int_] > -(dkwd("STEP")[int_ > int_]) ;
      rowsite_body = ctype > int_ > int_ > orient > -siterpt_stmt > ';' ;
      row = ctype > rowsite_body ;
      site = rowsite_body[at_c<1>(_val) = _1] ;

      dbu = dkwd("UNITS", 1)[dkwd("DISTANCE", 1)[dkwd("MICRONS", 1)[int_]]] > ';' ;

      semi_terminated = no_skip[*(char_ - ';') > ';'] ;
      def_file = dkwd("DESIGN")[dname[at_c<0>(_val) = _1]] > ';' >
	           // one giant case statement, of sorts, made possible by "dkwd"
	           (dkwd("VERSION", 0, 1)[(double_ > ';')[at_c<1>(_val) = _1]] /
                    dkwd("HISTORY")[semi_terminated [push_back(at_c<7>(_val), _1)]] /
		    dkwd("DIEAREA", 0, 1)[(rect[at_c<2>(_val) = _1] > ';')] /
		    dkwd("UNITS", 0, 1)[dkwd("DISTANCE")[dkwd("MICRONS")[(int_[at_c<3>(_val) = _1] > ';')]]] /
		    // a semi-manual approach to repeated components.  May need to make my own directive
		    // for cleanest approach.  You can pass in the rule as an inherited attribute,
		    // but the type of the parent rule needs to be different for each since the synthesized
		    // attribute is different (i.e., std::vector<attribute_of_child_rule>)
		    dkwd("COMPONENTS", 0, 1)[omit[int_[_a = _1]] > ';' >
					     repeat(_a)[component[push_back(at_c<4>(_val), _1)]] >
					     dkwd("END")[dkwd("COMPONENTS")[eps]]] /
		    dkwd("NETS", 0, 1)[omit[int_[_a = _1]] > ';' >
				       repeat(_a)[net[push_back(at_c<5>(_val), _1)]] >
				       dkwd("END")[dkwd("NETS")[eps]]] /
		    dkwd("ROW")[row[push_back(at_c<6>(_val), _1)]] /
		    dkwd("SITE")[site[push_back(at_c<6>(_val), _1)]] /
		    // catchall parsers for stuff we don't handle yet
		    dkwd((string("VIAS")|string("NETS")|string("SPECIALNETS")|string("PINS"))[_b=_1]
                       )[int_ > ';' > *('-' > semi_terminated) > dkwd("END", 1)[lit(_b)]] /
		    dkwd(string("TRACKS")|string("GCELLGRID"))[semi_terminated]) >
                 dkwd("END", 1)[kwd("DESIGN", 1)[eps]] ;

      // Debugging assistance

      dname.name("Design Name");
      version_stmt.name("VERSION");
      diearea_stmt.name("DIEAREA");
      component.name("Component");
      orient.name("Orientation");
      ctype.name("Cell Type");
      semi_terminated.name("Semicolon-terminated string");

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

  // helpful abbreviations
  typedef lefdefskipper<Iterator> skipper;

  template<typename Signature>
    struct Rule
    {
      typedef boost::spirit::qi::rule<Iterator, Signature, skipper> type;
    };

  // VERSION takes no parameters (a.k.a. "inherited attributes") and synthesizes a double for its attribute
  typename Rule<double()>::type version_stmt;

  // for HISTORY, and for grabbing things we don't yet parse
  typename Rule<std::string()>::type semi_terminated;

  // points "( x y )" produces defpoint structs (see deftypes.h)
  typename Rule<defpoint()>::type point;
  // rects "( llx lly ) ( urx ury )" synthesize defrect structs containing two points
  typename Rule<defrect()>::type rect, diearea_stmt;

  typedef typename Rule<std::string()>::type StringRule;
  StringRule orient, dname, ctype;

  // Site (or named row) statement
  typename Rule<siterepeat()>::type siterpt_stmt;
  typename Rule<rowsite()>::type row, site;
  typename Rule<rowsite_b()>::type rowsite_body;

  // storage for component names for ease of parse checking
  comp_symtab_t comp_symtab;

  // a single instance within the COMPONENTS section (name, celltype, placement)
  comp_parser<Iterator> component;

  // a single NET
  net_parser<Iterator> net;

  typename Rule<int()>::type dbu;

  // a catchall boost::spirit::qi::rule for everything I don't (yet) parse.  No attribute synthesized.
  boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<std::string>, skipper>
    tracks_stmt, gcellgrid_stmt, history_stmt;

  // The DEF file as a whole
  boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<int, std::string>,
    def(), skipper> def_file;

};

}

#endif
