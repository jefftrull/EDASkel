// A Boost Spirit-based LEF parser, part of EDASkel, a sample EDA app
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

#if !defined(EDASKEL_LEF_PARSER)
#define EDASKEL_LEF_PARSER

#include <vector>
#include <iostream>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/lexical_cast.hpp>

#include "leftypes.h"
#include "lefdef.h"



using namespace EDASkel;

namespace LefParse {

// just enough LEF to get library cells and their outlines

// Tokens
enum TokenIds {
  T_NAMESCASESENSITIVE = 1000, T_ON, T_OFF,
  T_LAYER, T_END, T_VIA, T_VIARULE, T_CLASS, T_PAD, T_CORE, T_ENDCAP,
  T_SITE, T_SYMMETRY, T_BY, T_SPACING, T_UNITS, T_FOREIGN, T_ORIGIN,
  T_SIZE, T_PIN, T_OBS, T_MACRO, T_LIBRARY,
  T_COVER, T_RING, T_BLOCK,
  T_INPUT, T_OUTPUT, T_SPACER, T_INOUT, T_POWER,
  T_FEEDTHRU, T_TIEHIGH, T_TIELOW,
  T_PRE, T_POST, T_TOPLEFT, T_TOPRIGHT, T_BOTTOMLEFT, T_BOTTOMRIGHT,
  T_X, T_Y, T_R90,
  T_ANY
};

template <typename Lexer>
struct LefTokens : boost::spirit::lex::lexer<Lexer>
{
  LefTokens()
    : double_("-?[0-9]+(\\.[0-9]+)?")
  {
    namespace lex = boost::spirit::lex;
    using boost::phoenix::ref;

    // Expecting simpler identifiers than in the DEF case
    ident = "[a-zA-Z][a-zA-Z0-9_]+";

    this->self =
      lex::string("NAMESCASESENSITIVE", T_NAMESCASESENSITIVE)
      | lex::string("ON", T_ON)
      | lex::string("OFF", T_OFF)
      | lex::string("LAYER", T_LAYER)
      | lex::string("END", T_END)
      | lex::string("VIA", T_VIA)
      | lex::string("VIARULE", T_VIARULE)
      | lex::string("CLASS", T_CLASS)
      | lex::string("PAD", T_PAD)
      | lex::string("CORE", T_CORE)
      | lex::string("ENDCAP", T_ENDCAP)
      | lex::string("SITE", T_SITE)
      | lex::string("SYMMETRY", T_SYMMETRY)
      | lex::string("BY", T_BY)
      | lex::string("SPACING", T_SPACING)
      | lex::string("UNITS", T_UNITS)
      | lex::string("FOREIGN", T_FOREIGN)
      | lex::string("ORIGIN", T_ORIGIN)
      | lex::string("SIZE", T_SIZE)
      | lex::string("PIN", T_PIN)
      | lex::string("OBS", T_OBS)
      | lex::string("MACRO", T_MACRO)
      | lex::string("LIBRARY", T_LIBRARY)
      | lex::string("COVER", T_COVER)
      | lex::string("RING", T_RING)
      | lex::string("BLOCK", T_BLOCK)
      | lex::string("INPUT", T_INPUT)
      | lex::string("OUTPUT", T_OUTPUT)
      | lex::string("SPACER", T_SPACER)
      | lex::string("INOUT", T_INOUT)
      | lex::string("POWER", T_POWER)
      | lex::string("FEEDTHRU", T_FEEDTHRU)
      | lex::string("TIEHIGH", T_TIEHIGH)
      | lex::string("TIELOW", T_TIELOW)
      | lex::string("PRE", T_PRE)
      | lex::string("POST", T_POST)
      | lex::string("TOPLEFT", T_TOPLEFT)
      | lex::string("TOPRIGHT", T_TOPRIGHT)
      | lex::string("BOTTOMLEFT", T_BOTTOMLEFT)
      | lex::string("BOTTOMRIGHT", T_BOTTOMRIGHT)
      | lex::string("X", T_X)
      | lex::string("Y", T_Y)
      | lex::string("R90", T_R90)
      | ident
      | double_
      | '+' | '-' | '(' | ')' | ';'
      ;

    // whitespace and comments
    this->self += lex::string("([ \\t\\n]+)|(#[^\\n]*\\n)")
        [
	 lex::_pass = lex::pass_flags::pass_ignore
        ];

    // catchall (mostly for error processing)
    this->self += lex::string(".", T_ANY);
  }

  // identifiers have a string attribute
  boost::spirit::lex::token_def<std::string> ident;
  // numbers
  boost::spirit::lex::token_def<double> double_;
 };

// for error handling
struct error_info_impl
{
  // required by phoenix::function; gives signature
  template <typename, typename>
  struct result
  {
    typedef std::string type;
  };

  template<typename Iterator, typename What>
  std::string operator()(Iterator const& actual_token_iter,
			 What const& what) const
  {
    return "Error! Expecting: " + boost::lexical_cast<std::string>(what) +
      " but saw: " + std::string(actual_token_iter->matched().begin(),
				 actual_token_iter->matched().end());
  }
};

template <typename Iterator, typename Lexer>
struct lefparser : boost::spirit::qi::grammar<Iterator, lef()>
{

  template <typename TokenDef>
  lefparser(TokenDef const& tok) : lefparser::base_type(lef_file)
    {
      using namespace boost::spirit::qi;
      namespace qi = boost::spirit::qi;         // BOZO silly to have both this and the previous
      using boost::spirit::_1;                  // access attributes for component count check
      using boost::phoenix::val;                // for error handling
      using boost::phoenix::construct;          // for error handling
      using boost::phoenix::at_c;               // to refer to pieces of wrapped structs
      using boost::phoenix::push_back;          // to store results in containers
      using boost::phoenix::ref;                // for storing parsing info into parser var (instead of attributes)

      // top-level elements in a LEF file

      casesens = qi::raw_token(T_NAMESCASESENSITIVE) >
	(qi::raw_token(T_ON)[ref(case_on) = true] |
	 qi::raw_token(T_OFF)[ref(case_on) = false]) > ';' ; 

      // general identifiers
      // celltypes: assuming only underscore might be used out of the non-alphanumeric chars
      // TBD with this lexer-based implementation, if we want to distinguish between acceptable
      // identifiers for e.g. nets and cells, we need to do check here somehow...
      // for now will just use the "ident" token

      layer = qi::raw_token(T_LAYER) > tok.ident[_a = _1] >> *!qi::raw_token(T_END) > qi::raw_token(T_END) > tok.ident[_b = _1] > eps(_a == _b) ;

      via = qi::raw_token(T_VIA) > tok.ident[_a = _1] >> *!qi::raw_token(T_END) >> qi::raw_token(T_END) > tok.ident[_b = _1] > eps(_a == _b) ;

      viarule = qi::raw_token(T_VIARULE) > tok.ident[_a = _1] >> *!qi::raw_token(T_END) >> qi::raw_token(T_END) > tok.ident[_b = _1] > eps(_a == _b) ;

      // BOZO review case sensitivity of these keywords and correct if necessary
      // BOZO this can appear in DEF as well; move to common LEF/DEF header
      // but how? this is internal to class.  Oh hey, inheritance...?

      // one-symbol site classes
      siteclass = (qi::raw_token(T_COVER) >> attr(SITECLASS_COVER)) |
	          (qi::raw_token(T_RING) >> attr(SITECLASS_RING)) |
	          (qi::raw_token(T_PAD) >> attr(SITECLASS_PAD)) |
	          (qi::raw_token(T_BLOCK) >> attr(SITECLASS_BLOCK)) |
	          (qi::raw_token(T_CORE) >> attr(SITECLASS_CORE));
      // classes distinguished by second symbol (all ENDCAP and some PAD, CORE)
      padclass =  (qi::raw_token(T_INPUT) >> attr(SITECLASS_PAD_INPUT)) |
	          (qi::raw_token(T_OUTPUT) >> attr(SITECLASS_PAD_OUTPUT)) |
	          (qi::raw_token(T_SPACER) >> attr(SITECLASS_PAD_SPACER)) |
	          (qi::raw_token(T_INOUT) >> attr(SITECLASS_PAD_INOUT)) |
	          (qi::raw_token(T_POWER) >> attr(SITECLASS_PAD_POWER));
      coreclass = (qi::raw_token(T_FEEDTHRU) >> attr(SITECLASS_CORE_FEEDTHRU)) |
	          (qi::raw_token(T_TIEHIGH) >> attr(SITECLASS_CORE_TIEHIGH)) |
	          (qi::raw_token(T_TIELOW) >> attr(SITECLASS_CORE_TIELOW));
      endcapclass = (qi::raw_token(T_PRE) >> attr(SITECLASS_ENDCAP_PRE)) ||
	          (qi::raw_token(T_POST) >> attr(SITECLASS_ENDCAP_POST)) ||
	          (qi::raw_token(T_TOPLEFT) >> attr(SITECLASS_ENDCAP_TOPLEFT)) ||
	          (qi::raw_token(T_TOPRIGHT) >> attr(SITECLASS_ENDCAP_TOPRIGHT)) ||
	          (qi::raw_token(T_BOTTOMLEFT) >> attr(SITECLASS_ENDCAP_BOTTOMLEFT)) ||
	          (qi::raw_token(T_BOTTOMRIGHT) >> attr(SITECLASS_ENDCAP_BOTTOMRIGHT));

      classrule = qi::raw_token(T_CLASS) >
	          ((qi::raw_token(T_PAD) >> padclass) |    // three kinds of two-keyword sites
		   (qi::raw_token(T_CORE) >> coreclass) |
		   (qi::raw_token(T_ENDCAP) > endcapclass) |
		   siteclass) > ';' ;               // try single keyword sites last

      sitesym = (qi::raw_token(T_X) >> attr(SITESYM_X)) |
	        (qi::raw_token(T_Y) >> attr(SITESYM_Y)) |
	        (qi::raw_token(T_R90) >> attr(SITESYM_R90));

      // BOZO add encountered sites to symbol table for checking
      // BOZO need a better way to handle case-insensitive keywords
      siterule %= qi::raw_token(T_SITE) > tok.ident[_a = _1] >
	classrule >
	-(qi::raw_token(T_SYMMETRY) > +sitesym > ';' ) >
	// BOZO use lefextent here
	qi::raw_token(T_SIZE) > tok.double_ > qi::raw_token(T_BY) > tok.double_ > ';' >
	// BOZO optional ROWPATTERN goes here
	qi::raw_token(T_END) > omit[tok.ident[_b = _1]] > eps[_a == _b];

      spacing = qi::raw_token(T_SPACING) > *!qi::raw_token(T_END) >> qi::raw_token(T_END) > qi::raw_token(T_SPACING) ;
      units = qi::raw_token(T_UNITS) > *!qi::raw_token(T_END) > qi::raw_token(T_END) >> qi::raw_token(T_UNITS) ;

      // macro properties
      point = tok.double_ >> tok.double_ ;
      foreign = qi::raw_token(T_FOREIGN) > tok.ident > point > ';' ;
      origin  = qi::raw_token(T_ORIGIN) > point > ';' ;
      macrosymmetry = qi::raw_token(T_SYMMETRY) > *sitesym > ';' ;
      macrosize = qi::raw_token(T_SIZE) > tok.double_ > no_case[qi::raw_token(T_BY)] > tok.double_ > ';' ;

      // macro elements
      pin = qi::raw_token(T_PIN) > tok.ident[_a = _1] > *!qi::raw_token(T_END) > qi::raw_token(T_END) > tok.ident[_b = _1] > eps(_a == _b) ;
      obs = qi::raw_token(T_OBS) > *!qi::raw_token(T_END) > qi::raw_token(T_END) ;

      macro = qi::raw_token(T_MACRO) > tok.ident[_a = _1, at_c<0>(_val) = _1] >
	// these can be in any order, but the caret checks that at least one is present
	// SIZE is required, so I guess that makes it work
	// have to use semantic actions here because the () grouped stuff becomes its own sequence
	(classrule[at_c<1>(_val) = _1] ^
	 foreign[at_c<2>(_val) = _1] ^
	 origin[at_c<3>(_val) = _1] ^
	 macrosize[at_c<4>(_val) = _1] ^
	 macrosymmetry[at_c<5>(_val) = _1] ^
	 (qi::raw_token(T_SITE) > tok.ident[at_c<6>(_val) = _1] > ';')) >
	// BOZO not stored for now:
	*pin >
	-obs >
	qi::raw_token(T_END) > tok.ident[_b = _1] > eps(_a == _b) ;


      // define some major elements

      // strangely a LEF file does not begin with anything distinctive... just goes right into the data
      // BOZO siterule should store site name into symtab
      // evidently EVERY statement is optional
      lef_file = *(casesens | units | spacing | layer | via | viarule |
		   siterule [push_back(at_c<0>(_val), _1)] |
		   macro [push_back(at_c<1>(_val), _1)])
	>> -(qi::raw_token(T_END) > qi::raw_token(T_LIBRARY)) ;

      // Debugging assistance
      BOOST_SPIRIT_DEBUG_NODE(casesens);
      BOOST_SPIRIT_DEBUG_NODE(layer);
      BOOST_SPIRIT_DEBUG_NODE(via);
      BOOST_SPIRIT_DEBUG_NODE(viarule);
      BOOST_SPIRIT_DEBUG_NODE(siterule);
      BOOST_SPIRIT_DEBUG_NODE(spacing);
      BOOST_SPIRIT_DEBUG_NODE(units);
      BOOST_SPIRIT_DEBUG_NODE(macro);
      BOOST_SPIRIT_DEBUG_NODE(foreign);
      BOOST_SPIRIT_DEBUG_NODE(origin);
      BOOST_SPIRIT_DEBUG_NODE(point);
      BOOST_SPIRIT_DEBUG_NODE(classrule);
      BOOST_SPIRIT_DEBUG_NODE(padclass);
      BOOST_SPIRIT_DEBUG_NODE(coreclass);
      BOOST_SPIRIT_DEBUG_NODE(endcapclass);
      BOOST_SPIRIT_DEBUG_NODE(lef_file);
      BOOST_SPIRIT_DEBUG_NODE(pin);
      BOOST_SPIRIT_DEBUG_NODE(obs);

      on_error<fail>
        (
	 lef_file
	 , std::cerr << error_info(boost::spirit::_3, boost::spirit::_4) << std::endl
	 );
    }

  // define rules (inherited attributes (a.k.a. parameters), local variables, synthesized attributes (a.k.a. return values )
  // globals for the parser object
  bool case_on;   // case sensitivity ON/OFF -> true/false

  // rules neither inheriting nor synthesizing an attribute
  typedef boost::spirit::qi::rule<Iterator> NoAttrRule;
  NoAttrRule spacing, units, casesens, obs;

  // for now macro synthesizes a string AND has a local variable (for checking END statement)
  boost::spirit::qi::rule<Iterator, lefmacro(),
			  boost::spirit::qi::locals<std::string, std::string> > macro;

  // Some keywords turn directly into enums in the database
  boost::spirit::qi::rule<Iterator, SiteClass()> siteclass, padclass, coreclass, endcapclass;
  boost::spirit::qi::rule<Iterator, SiteSymmetry()> sitesym;
  boost::spirit::qi::rule<Iterator, SiteClass()> classrule;
  // site rule uses local var (to check END statement) and synthesizes a Site
  boost::spirit::qi::rule<Iterator, Site(),
			  boost::spirit::qi::locals<std::string, std::string> > siterule;
  // macros can list their site symmetry options
  boost::spirit::qi::rule<Iterator, std::vector<SiteSymmetry>()> macrosymmetry;
  // macros give "origin" as a real numbered point
  boost::spirit::qi::rule<Iterator, lefpoint()> point, origin;
  // size (width/height) is identical in terms of data stored but has a different interpretation
  boost::spirit::qi::rule<Iterator, lefextent()> macrosize;
  // macro FOREIGN statement gives a name and a point.  BOZO is it legal to name something other than your own macro?
  boost::spirit::qi::rule<Iterator, lefforeign()> foreign;

  // PIN and OBS are not currently handled.  OBS will be necessary for polygonal block outlines. PIN needed when we do nets

  // a catchall rule for everything I don't (yet) parse.  No attribute synthesized, but it has a local string var
  // (for checking END statements)
  typedef boost::spirit::qi::rule<Iterator,
				  boost::spirit::qi::locals<std::string, std::string> > Unparsed;
  Unparsed layer, via, viarule, pin;

  // stuff inside macros I don't handle
  boost::spirit::qi::rule<Iterator > unparsed_macro_stuff;

  // The LEF file as a whole
  boost::spirit::qi::rule<Iterator, lef()> lef_file;

  // error handler
  boost::phoenix::function<error_info_impl> error_info;
};

}

#endif
