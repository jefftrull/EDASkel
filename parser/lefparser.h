// A Boost Spirit-based LEF parser, part of EDASkel, a sample EDA app
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

#if !defined(EDASKEL_LEF_PARSER)
#define EDASKEL_LEF_PARSER

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/home/phoenix/statement/sequence.hpp>

#include "keyword.h"
#include "leftypes.h"

#include <vector>
#include <iostream>

namespace LefParse {

// A skip parser for LEF comments and spaces
// adapted from a presentation at Boostcon 2010 by Michael Caisse
template <typename Iterator>
struct lefskipper : boost::spirit::qi::grammar< Iterator >
{
 lefskipper() : lefskipper::base_type(skip_it)
    {
      using namespace boost::spirit::qi;

      comment = '#' >> *( char_ - eol ) >> eol ;
      skip_it = comment | space ;
    }
  boost::spirit::qi::rule<Iterator> skip_it;
  boost::spirit::qi::rule<Iterator> comment;
};


// just enough LEF to get library cells and their outlines
template <typename Iterator>
struct lefparser : boost::spirit::qi::grammar<Iterator,
                                              lef(),
                                              lefskipper<Iterator> >
{

  lefparser() : lefparser::base_type(lef_file)
    {
      using namespace boost::spirit::qi;
      using namespace distinct;                 // for keywords.  We are a "phrase" (not character) parser
                                                // so need to distinguish keywords from following alphanumerics
      using boost::spirit::_1;                  // access attributes for component count check
      using boost::phoenix::val;                // for error handling
      using boost::phoenix::construct;          // for error handling
      using boost::phoenix::at_c;               // to refer to pieces of wrapped structs
      using boost::phoenix::ref;                // for storing parsing info into parser var (instead of attributes)

      // top-level elements in a LEF file

      casesens = keyword["NAMESCASESENSITIVE"] >
	         omit[string("ON")[ref(case_on) = true] |
		      string("OFF")[ref(case_on) = false]] > ';' ; 

      // general identifiers
      // celltypes: assuming only underscore might be used out of the non-alphanumeric chars
      ctype %= lexeme[alpha >> *(alnum | char_('_'))] ;
      id %= lexeme[alpha >> *(alnum | char_('-') | char_('_'))] ;

      layer = keyword["LAYER"] > id[_a = _1] > *(char_ - (keyword["END"] >> lit(_a))) > keyword["END"] > lit(_a) ;

      via = keyword["VIA"] > id[_a = _1] > *(char_ - (keyword["END"] >> lit(_a))) > keyword["END"] > lit(_a) ;

      viarule = keyword["VIARULE"] > id[_a = _1] > *(char_ - (keyword["END"] >> lit(_a))) > keyword["END"] > lit(_a) ;

      // BOZO review case sensitivity of these keywords and correct if necessary
      // BOZO this can appear in DEF as well; move to common LEF/DEF header
      // but how? this is internal to class.  Oh hey, inheritance...?

      // one-symbol site classes
      siteclass.add("cover", SITECLASS_COVER)("ring", SITECLASS_RING)("pad", SITECLASS_PAD);
      siteclass.add("block", SITECLASS_BLOCK)("core", SITECLASS_CORE);
      // classes distinguished by second symbol (all ENDCAP and some PAD, CORE)
      padclass.add("input", SITECLASS_PAD_INPUT)("output", SITECLASS_PAD_OUTPUT)("spacer", SITECLASS_PAD_SPACER);
      padclass.add("inout", SITECLASS_PAD_INOUT)("power", SITECLASS_PAD_POWER);
      coreclass.add("feedthru", SITECLASS_CORE_FEEDTHRU)("tiehigh", SITECLASS_CORE_TIEHIGH);
      coreclass.add("tielow", SITECLASS_CORE_TIELOW);
      endcapclass.add("pre", SITECLASS_ENDCAP_PRE)("post", SITECLASS_ENDCAP_POST)("topleft", SITECLASS_ENDCAP_TOPLEFT);
      endcapclass.add("topright", SITECLASS_ENDCAP_TOPRIGHT)("bottomleft", SITECLASS_ENDCAP_BOTTOMLEFT)("bottomright", SITECLASS_ENDCAP_BOTTOMRIGHT);

      // note there may be some more clever 2-level "Nabialek" style approach that will work
      classrule = keyword["CLASS"] >
	((no_case[keyword["pad"]] >> no_case[keyword[padclass]]) |   // three kinds of two-keyword sites
	 (no_case[keyword["core"]] >> no_case[keyword[coreclass]]) |
	 (no_case[keyword["endcap"]] >> no_case[keyword[endcapclass]]) |
	 no_case[keyword[siteclass]]) > ';' ;               // try single keyword sites last

      sitesym.add("x", SITESYM_X)("y", SITESYM_Y)("r90", SITESYM_R90);
      // BOZO add encountered sites to symbol table for checking
      // BOZO need a better way to handle case-insensitive keywords
      siterule %= (keyword["SITE"] | keyword["site"] ) > id[_a = _1] >
	classrule >
	-(keyword["SYMMETRY"] > +no_case[keyword[sitesym]] > ';' ) >
	// BOZO use lefextent here
	no_case[keyword["size"]] > float_ > no_case[keyword["by"]] > float_ > ';' >
	// BOZO optional ROWPATTERN goes here
	no_case[keyword["end"]] > omit[string(_a)] ;

      spacing = keyword["SPACING"] > *(char_ - (keyword["END"] >> keyword["SPACING"])) > keyword["END"] > keyword["SPACING"] ;
      units = keyword["UNITS"] > *char_ > keyword["END"] > keyword["UNITS"] ;

      // macro properties
      point = float_ >> float_ ;
      foreign = keyword["FOREIGN"] > ctype > point > ';' ;
      origin  = keyword["ORIGIN"] > point > ';' ;
      macrosymmetry = keyword["SYMMETRY"] > *no_case[keyword[sitesym]] > ';' ;
      macrosize = keyword["SIZE"] > float_ > no_case[keyword["by"]] > float_ > ';' ;

      // macro elements
      pin = keyword["PIN"] > id[_a = _1] > omit[*(char_ - (keyword["END"] >> keyword[lit(_a)]))] > keyword["END"] > keyword[lit(_a)] ;
      obs = keyword["OBS"] > omit[*(!keyword["END"] >> char_)] > keyword["END"] ;

      macro = keyword["MACRO"] > ctype[_a = _1, at_c<0>(_val) = _1] >
	// these can be in any order, but the caret checks that at least one is present
	// SIZE is required, so I guess that makes it work
	// have to use semantic actions here because the () grouped stuff becomes its own sequence
	(classrule[at_c<1>(_val) = _1] ^
	 foreign[at_c<2>(_val) = _1] ^
	 origin[at_c<3>(_val) = _1] ^
	 macrosize[at_c<4>(_val) = _1] ^
	 macrosymmetry[at_c<5>(_val) = _1] ^
	 (keyword["SITE"] > id[at_c<6>(_val) = _1] > ';')) >
	// BOZO not stored for now:
	*pin >
	-obs >
	keyword["END"] > lit(_a) ;


      // define some major elements

      // strangely a LEF file does not begin with anything distinctive... just goes right into the data
      // BOZO verify that ^ means optional and if not, do appropriate coding
      // BOZO siterule should store site name into symtab
      //      lef_file %= (casesens ^ units ^ spacing ^ *layer ^ *via ^ *viarule ^ *siterule ^ *macro ) || (keyword["END"] > keyword["LIBRARY"]) ;
      lef_file %= -casesens >> -units >> *layer >> *via >> *viarule >> -spacing >> *siterule >> *macro >> -(keyword["END"] > keyword["LIBRARY"]) ;

      // Debugging assistance
      casesens.name("NAMESCASESENSITIVE");
      ctype.name("Cell type");
      id.name("Identifier");
      layer.name("Layer definition");
      via.name("Via definition");
      viarule.name("Via Rule");
      siterule.name("Site definition");
      spacing.name("Global spacing definition");
      units.name("Units definition");
      macro.name("Macro definition");
      lef_file.name("Top-level LEF file");
      pin.name("PIN definition");
      obs.name("OBS definition");

      on_error<fail>
        (
	 lef_file
	 , std::cerr
	 << val("Error! Expecting ")
	 << boost::spirit::_4                               // what failed?
	 << val(" here: \"")
	 << construct<std::string>(boost::spirit::_3, boost::spirit::_2)   // iterators to error-pos, end
	 << val("\"")
	 << std::endl
	 );
    }

  // define rules (inherited attributes (a.k.a. parameters), local variables, synthesized attributes (a.k.a. return values )
  // globals for the parser object
  bool case_on;   // case sensitivity ON/OFF -> true/false

  // rules neither inheriting nor synthesizing an attribute
  typedef boost::spirit::qi::rule<Iterator, lefskipper<Iterator> > NoAttrRule;
  NoAttrRule spacing, units, casesens, obs;

  // rules synthesizing a string
  typedef boost::spirit::qi::rule<Iterator, std::string(), lefskipper<Iterator> > StringRule;
  StringRule ctype, id;

  // for now macro synthesizes a string AND has a local variable (for checking END statement)
  boost::spirit::qi::rule<Iterator, lefmacro(), boost::spirit::qi::locals<std::string>, lefskipper<Iterator> > macro;

  // Some keywords turn directly into enums in the database
  boost::spirit::qi::symbols<char, SiteClass> siteclass, padclass, coreclass, endcapclass;
  boost::spirit::qi::symbols<char, SiteSymmetry> sitesym;
  boost::spirit::qi::rule<Iterator, SiteClass(), lefskipper<Iterator> > classrule;
  // site rule uses local var (to check END statement) and synthesizes a Site
  boost::spirit::qi::rule<Iterator, Site(), boost::spirit::qi::locals<std::string>, lefskipper<Iterator>  > siterule;
  // macros can list their site symmetry options
  boost::spirit::qi::rule<Iterator, std::vector<SiteSymmetry>(), lefskipper<Iterator> > macrosymmetry;
  // macros give "origin" as a real numbered point
  boost::spirit::qi::rule<Iterator, lefpoint(), lefskipper<Iterator> > point, origin;
  // size (width/height) is identical in terms of data stored but has a different interpretation
  boost::spirit::qi::rule<Iterator, lefextent(), lefskipper<Iterator> > macrosize;
  // macro FOREIGN statement gives a name and a point.  BOZO is it legal to name something other than your own macro?
  boost::spirit::qi::rule<Iterator, lefforeign(), lefskipper<Iterator> > foreign;

  // PIN and OBS are not currently handled.  OBS will be necessary for polygonal block outlines. PIN needed when we do nets

  // a catchall rule for everything I don't (yet) parse.  No attribute synthesized, but it has a local string var
  // (for checking END statements)
  typedef boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<std::string>, lefskipper<Iterator>  > Unparsed;
  Unparsed layer, via, viarule, pin;

  // stuff inside macros I don't handle
  boost::spirit::qi::rule<Iterator, lefskipper<Iterator>  > unparsed_macro_stuff;

  // The LEF file as a whole
  boost::spirit::qi::rule<Iterator, lef(), lefskipper<Iterator> > lef_file;


};

}

#endif
