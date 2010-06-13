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
#if !defined(EDASKEL_DEF_PARSER)
#define EDASKEL_DEF_PARSER

#include <vector>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include "keyword.h"
#include "deftypes.h"
#include "lefdef.h"

using namespace EDASkel;

namespace DefParse {

// a starter DEF grammar
template <typename Iterator>
struct defparser : boost::spirit::qi::grammar<Iterator,
                                              def(),
                                              lefdefskipper<Iterator> >
{

  defparser() : defparser::base_type(def_file)
    {
      using namespace boost::spirit::qi;
      using namespace distinct;                 // for keywords.  We are a "phrase" (not character) parser
                                                // so need to distinguish keywords from following alphanumerics
      using boost::spirit::_1;                  // access attributes for component count check
      using boost::phoenix::push_back;          // to store results in containers
      using boost::phoenix::val;                // for error handling
      using boost::phoenix::construct;          // for error handling
      using boost::phoenix::at_c;               // to refer to pieces of wrapped structs

      // top-level elements in a DEF file
      version_stmt = keyword["VERSION"] > double_ > ';' ;

      point %= '(' >> int_ >> int_ >> ')' ;       // points are parenthesized pairs, no comma
      rect %= point >> point ;                    // rects are just two points in a row
      diearea_stmt %= keyword["DIEAREA"] > rect > ';' ;

      // identifiers
      // design name: letter followed by letters, numbers, underscore, or hyphen
      dname %= lexeme[alpha >> *(alnum | char_('-') | char_('_'))] ;
      // instance names: letter followed by letters, numbers, underscore, hyphen, square brackets, slashes (hierarchy)
      // BOZO support escapes
      iname %= lexeme[alpha >> *(alnum | char_('-') | char_('_') | char_('[') | char_(']') | char_('/'))] ;
      // celltypes: assuming only underscore might be used out of the non-alphanumeric chars
      ctype %= lexeme[alpha >> *(alnum | char_('_'))] ;

      // define some major elements

      // components (instances)
      orient = keyword[repeat(1)[char_("NSEW")] |         // a single cardinal direction, OR
		       (char_('F') >> char_("NSEW"))] ;   // the same, but flipped

      // Using ">" here instead of ">>" implies a required sequence and allows the parser to check it
      // specifically (instead of simply failing on the whole input)
      plcinfo %= '+' >> (keyword[string("FIXED")] | keyword[string("PLACED")]) >
		     point > orient ;    // location and orientation

      weight %= '+' >> keyword["WEIGHT"] > int_ ;

      source = '+' >> keyword["SOURCE"] > (keyword["DIST"] | keyword["NETLIST"] | keyword["USER"] | keyword["TIMING"]) ;

      // components required instance name and celltype; optional any (or none) of placement and weight, in any order:
      component %= '-' > iname > ctype > (plcinfo ^ omit[weight] ^ source ^ eps ) > ';' ;

      // define repeat rules for each element section
      // I tried to make this generic but failed.  You can pass in the rule as an inherited attribute,
      // but the type of the parent rule needs to be different for each since the synthesized attribute is
      // different (i.e., std::vector<attribute_of_child_rule>)
      comps_section %= keyword["COMPONENTS"] > omit[int_[_a = _1]] > ';' >  // remember count in a local variable
	               repeat(_a)[component] >                        // expect that many copies of the supplied rule
	               keyword["END"] > keyword["COMPONENTS"] ;       // END followed by the section name again

      // My copy of the LEF/DEF reference does not show this SITE command as valid for DEF yet my example data does...
      // The example data's syntax is very similar to that defined for ROW, so I'll combine them
      siterpt_stmt = keyword["DO"] > int_ > keyword["BY"] > int_ > -(keyword["STEP"] > int_ > int_) ;
      rowsite_stmt = ((keyword["ROW"] > ctype) | keyword["SITE"]) > ctype > int_ > int_ > orient >
	             -siterpt_stmt > ';' ;

      dbu = keyword["UNITS"] > keyword["DISTANCE"] > keyword["MICRONS"] > int_ > ';' ;

      // This parser only handles components and a couple of misc. statements
      // here's a catchall parser to discard all other data
      // BOZO update this for all keywords we parse - make a symbol table out of them for ease of use and speed
      tracks_stmt = keyword["TRACKS"] > *(char_ - ';') > ';' ;
      gcellgrid_stmt = keyword["GCELLGRID"] > *(char_ - ';') > ';' ;
      history_stmt = keyword["HISTORY"] > *(char_ - ';') > ';' ;

      unparsed = ((keyword[string("VIAS")] | keyword[string("NETS")] |
		   keyword[string("SPECIALNETS")] | keyword[string("PINS")])[_a = _1] > int_ > ';' >
		  *('-' > *(char_ - ';') > ';') > keyword["END"] > keyword[_a] ) | tracks_stmt | gcellgrid_stmt | history_stmt ;

      def_file = keyword["DESIGN"] > dname[at_c<0>(_val) = _1] > ';' >
                 *(version_stmt[at_c<1>(_val) = _1] |
		   diearea_stmt[at_c<2>(_val) = _1] |
		   dbu[at_c<3>(_val) = _1] |
      	           comps_section[at_c<4>(_val) = _1] |
		   rowsite_stmt[push_back(at_c<5>(_val), _1)] |
		   unparsed) >
      	         keyword["END"] > keyword["DESIGN"] ;

      // Debugging assistance

      dname.name("Design Name");
      iname.name("Instance Name");
      ctype.name("Cell Type");
      version_stmt.name("VERSION");
      diearea_stmt.name("DIEAREA");
      comps_section.name("COMPONENTS Section");
      component.name("Component");
      orient.name("Orientation");
      plcinfo.name("Optional Placement Info");


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

  // VERSION takes no parameters (a.k.a. "inherited attributes") and synthesizes a double for its attribute
  boost::spirit::qi::rule<Iterator, double(), lefdefskipper<Iterator> > version_stmt;

  // points "( x y )" produces defpoint structs (see deftypes.h)
  boost::spirit::qi::rule<Iterator, defpoint(), lefdefskipper<Iterator> > point;
  // rects "( llx lly ) ( urx ury )" synthesize defrect structs containing two points
  boost::spirit::qi::rule<Iterator, defrect(), lefdefskipper<Iterator> > rect, diearea_stmt;

  typedef boost::spirit::qi::rule<Iterator, std::string(), lefdefskipper<Iterator> > StringRule;
  StringRule orient, dname, iname, ctype;

  // rules neither inheriting nor synthesizing an attribute
  typedef boost::spirit::qi::rule<Iterator, lefdefskipper<Iterator> > NoAttrRule;
  NoAttrRule weight;

  // Site (or named row) statement
  boost::spirit::qi::rule<Iterator, siterepeat(), lefdefskipper<Iterator> > siterpt_stmt;
  boost::spirit::qi::rule<Iterator, rowsite(), lefdefskipper<Iterator> > rowsite_stmt;

  // rules pertaining to COMPONENTS - see deftypes.h for data results

  // optional placement info (placed vs. fixed, location, orientation)
  boost::spirit::qi::rule<Iterator, defplcinfo(), lefdefskipper<Iterator> > plcinfo;
  // a single instance within the COMPONENTS section (name, celltype, placement)
  boost::spirit::qi::rule<Iterator, defcomponent(), lefdefskipper<Iterator> > component;
  // a rule representing the entire COMPONENTS section
  boost::spirit::qi::rule<Iterator, std::vector<defcomponent>(), boost::spirit::qi::locals<int>, lefdefskipper<Iterator> > comps_section;

  boost::spirit::qi::rule<Iterator, int(), lefdefskipper<Iterator> > dbu;

  // a catchall rule for everything I don't (yet) parse.  No attribute synthesized.
  boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<std::string>, lefdefskipper<Iterator> > unparsed, tracks_stmt, gcellgrid_stmt, history_stmt, source;

  // The DEF file as a whole
  boost::spirit::qi::rule<Iterator, def(), lefdefskipper<Iterator> > def_file;

};

}

#endif
