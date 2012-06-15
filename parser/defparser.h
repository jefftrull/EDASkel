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
template <typename Iterator>
struct defparser : boost::spirit::qi::grammar<Iterator,
                                              boost::spirit::qi::locals<int, std::string>,
                                              def(),
                                              lefdefskipper<Iterator> >
{

  defparser() : defparser::base_type(def_file)
    {
      using namespace boost::spirit::qi;

      using boost::spirit::_1;                  // access attributes for component count check
      using boost::phoenix::push_back;          // to store results in containers
      using boost::phoenix::val;                // for error handling
      using boost::phoenix::construct;          // for error handling
      using boost::phoenix::at_c;               // to refer to pieces of wrapped structs
      using boost::spirit::repository::dkwd;
      using boost::spirit::repository::kwd;

      // top-level elements in a DEF file
      version_stmt = dkwd("VERSION", 1)[double_] > ';' ;

      point %= '(' >> int_ >> int_ >> ')' ;       // points are parenthesized pairs, no comma
      rect %= point >> point ;                    // rects are just two points in a row
      diearea_stmt %= dkwd("DIEAREA", 1)[rect] > ';' ;

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
      orient = repeat(1)[char_("NSEW")] |         // a single cardinal direction, OR
	        (char_('F') >> char_("NSEW")) ;   // the same, but flipped

      // Using ">" here instead of ">>" implies a required sequence and allows the parser to check it
      // specifically (instead of simply failing on the whole input)
      plcinfo %= '+' >> (dkwd("FIXED", 1)[attr("FIXED") > point > orient] |
			 dkwd("PLACED", 1)[attr("PLACED") > point > orient]) ;    // location and orientation

      weight %= '+' >> dkwd("WEIGHT", 1)[int_] ;

      source = '+' >> dkwd("SOURCE", 1) [(dkwd("DIST")[attr("DIST")] /
					  dkwd("NETLIST")[attr("NETLIST")] /
					  dkwd("USER")[attr("USER")] /
					  dkwd("TIMING")[attr("TIMING")])] ;

      // components required instance name and celltype; optional any (or none) of placement and weight, in any order:
      component %= '-' > iname > ctype > (plcinfo ^ omit[weight] ^ source ^ eps ) > ';' ;

      // My copy of the LEF/DEF reference does not show this SITE command as valid for DEF yet my example data does...
      // The example data's syntax is very similar to that defined for ROW, so I'll combine them
      siterpt_stmt = dkwd("DO")[int_] > dkwd("BY")[int_] > -(dkwd("STEP")[int_ > int_]) ;
      rowsite_body = ctype > int_ > int_ > orient > -siterpt_stmt > ';' ;
      row = ctype > rowsite_body ;
      site = rowsite_body[at_c<1>(_val) = _1] ;

      dbu = dkwd("UNITS", 1)[dkwd("DISTANCE", 1)[dkwd("MICRONS", 1)[int_]]] > ';' ;

      def_file = dkwd("DESIGN")[dname[at_c<0>(_val) = _1]] > ';' >
	           // one giant case statement, of sorts, made possible by "dkwd"
	           (dkwd("VERSION", 0, 1)[(double_ > ';')[at_c<1>(_val) = _1]] /
		    dkwd("DIEAREA", 0, 1)[(rect[at_c<2>(_val) = _1] > ';')] /
		    dkwd("UNITS", 0, 1)[dkwd("DISTANCE")[dkwd("MICRONS")[(int_[at_c<3>(_val) = _1] > ';')]]] /
		    // a semi-manual approach to repeated components.  May need to make my own directive
		    // for cleanest approach.  You can pass in the rule as an inherited attribute,
		    // but the type of the parent rule needs to be different for each since the synthesized
		    // attribute is different (i.e., std::vector<attribute_of_child_rule>)
		    dkwd("COMPONENTS", 0, 1)[omit[int_[_a = _1]] > ';' >
					     repeat(_a)[component[push_back(at_c<4>(_val), _1)]] >
					     dkwd("END")[dkwd("COMPONENTS")[eps]]] /
		    dkwd("ROW")[row[push_back(at_c<5>(_val), _1)]] /
		    dkwd("SITE")[site[push_back(at_c<5>(_val), _1)]] /
		    // catchall parsers for stuff we don't handle yet
		    dkwd((string("VIAS")|string("NETS")|string("SPECIALNETS")|string("PINS"))[_b=_1]
			 )[int_ > ';' > *('-' > *(char_ - ';') > ';') > dkwd("END", 1)[lit(_b)]] /
		    dkwd(string("TRACKS")|string("GCELLGRID")|string("HISTORY"))[*(char_ - ';') > ';']) >
	dkwd("END", 1)[kwd("DESIGN", 1)[eps]] ;

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
  boost::spirit::qi::rule<Iterator, rowsite(), lefdefskipper<Iterator> > row, site;
  boost::spirit::qi::rule<Iterator, rowsite_b(), lefdefskipper<Iterator> > rowsite_body;

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
  boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<int, std::string>,
    def(), lefdefskipper<Iterator> > def_file;

};

}

#endif
