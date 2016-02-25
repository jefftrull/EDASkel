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

#include <boost/spirit/include/lex_lexertl.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/lexical_cast.hpp>

#include "leftypes.h"


namespace LefParse {

// just enough LEF to get library cells and their outlines

// Tokens
enum TokenIds {
  T_NAMESCASESENSITIVE = 1000, T_ON, T_OFF, T_VERSION,
  T_LAYER, T_END, T_VIA, T_VIARULE, T_CLASS, T_PAD, T_CORE, T_ENDCAP,
  T_SITE, T_SYMMETRY, T_BY, T_SPACING, T_UNITS, T_FOREIGN, T_ORIGIN,
  T_SIZE, T_PIN, T_OBS, T_MACRO, T_LIBRARY, T_PORT,
  T_COVER, T_RING, T_BLOCK,
  T_INPUT, T_OUTPUT, T_SPACER, T_INOUT, T_POWER,
  T_FEEDTHRU, T_TIEHIGH, T_TIELOW,
  T_PRE, T_POST, T_TOPLEFT, T_TOPRIGHT, T_BOTTOMLEFT, T_BOTTOMRIGHT,
  T_X, T_Y, T_R90,
  T_TYPE, T_CUT, T_MASTERSLICE, T_IMPLANT, T_ROUTING, T_OVERLAP,
  T_ANY
};

template <typename Lexer>
struct LefTokens : boost::spirit::lex::lexer<Lexer>
{
  LefTokens()
    : double_("-?([0-9]+(\\.[0-9]+)?)|(\\.[0-9]+)")
  {
    namespace lex = boost::spirit::lex;
    using boost::phoenix::ref;

    // Expecting simpler identifiers than in the DEF case
    ident = "[a-zA-Z][a-zA-Z0-9_]*";

    this->self =
      lex::string("NAMESCASESENSITIVE", T_NAMESCASESENSITIVE)
      | lex::string("(?i:on)", T_ON)
      | lex::string("(?i:off)", T_OFF)
      | lex::string("(?i:version)", T_VERSION)
      | lex::string("(?i:layer)", T_LAYER)
      | lex::string("(?i:end)", T_END)
      | lex::string("(?i:via)", T_VIA)
      | lex::string("(?i:viarule)", T_VIARULE)
      | lex::string("(?i:class)", T_CLASS)
      | lex::string("(?i:pad)", T_PAD)
      | lex::string("(?i:core)", T_CORE)
      | lex::string("(?i:endcap)", T_ENDCAP)
      | lex::string("(?i:site)", T_SITE)
      | lex::string("(?i:symmetry)", T_SYMMETRY)
      | lex::string("(?i:by)", T_BY)
      | lex::string("(?i:spacing)", T_SPACING)
      | lex::string("(?i:units)", T_UNITS)
      | lex::string("(?i:foreign)", T_FOREIGN)
      | lex::string("(?i:origin)", T_ORIGIN)
      | lex::string("(?i:size)", T_SIZE)
      | lex::string("(?i:pin)", T_PIN)
      | lex::string("(?i:obs)", T_OBS)
      | lex::string("(?i:macro)", T_MACRO)
      | lex::string("(?i:library)", T_LIBRARY)
      | lex::string("(?i:port)", T_PORT)
      | lex::string("(?i:cover)", T_COVER)
      | lex::string("(?i:ring)", T_RING)
      | lex::string("(?i:block)", T_BLOCK)
      | lex::string("(?i:input)", T_INPUT)
      | lex::string("(?i:output)", T_OUTPUT)
      | lex::string("(?i:spacer)", T_SPACER)
      | lex::string("(?i:inout)", T_INOUT)
      | lex::string("(?i:power)", T_POWER)
      | lex::string("(?i:feedthru)", T_FEEDTHRU)
      | lex::string("(?i:tiehigh)", T_TIEHIGH)
      | lex::string("(?i:tielow)", T_TIELOW)
      | lex::string("(?i:pre)", T_PRE)
      | lex::string("(?i:post)", T_POST)
      | lex::string("(?i:topleft)", T_TOPLEFT)
      | lex::string("(?i:topright)", T_TOPRIGHT)
      | lex::string("(?i:bottomleft)", T_BOTTOMLEFT)
      | lex::string("(?i:bottomright)", T_BOTTOMRIGHT)
      | lex::string("(?i:x)", T_X)
      | lex::string("(?i:y)", T_Y)
      | lex::string("(?i:r90)", T_R90)
      | lex::string("(?i:type)", T_TYPE)
      | lex::string("(?i:cut)", T_CUT)
      | lex::string("(?i:masterslice)", T_MASTERSLICE)
      | lex::string("(?i:routing)", T_ROUTING)
      | lex::string("(?i:overlap)", T_OVERLAP)
      | lex::string("(?i:implant)", T_IMPLANT)
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
  template <typename Signature>
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
      using boost::phoenix::at_c;               // to refer to pieces of wrapped structs
      using boost::phoenix::push_back;          // to store results in containers
      using boost::phoenix::ref;                // for storing parsing info into parser var (instead of attributes)
      using qi::raw_token;                      // to turn token IDs into qi parsers

      // top-level elements in a LEF file

      casesens = raw_token(T_NAMESCASESENSITIVE) >
	(raw_token(T_ON)[ref(case_on) = true] |
	 raw_token(T_OFF)[ref(case_on) = false]) > ';' ; 

      version_stmt = raw_token(T_VERSION) > tok.double_ > ';' ;

      // general identifiers
      // celltypes: assuming only underscore might be used out of the non-alphanumeric chars
      // TBD with this lexer-based implementation, if we want to distinguish between acceptable
      // identifiers for e.g. nets and cells, we need to do check here somehow...
      // for now will just use the "ident" token

      // can't seem to successfully say "any token other than END" so I've got to do this:
      catchall = ';' | tok.ident | tok.double_ | raw_token(T_SPACING);

      layer = raw_token(T_LAYER) > tok.ident[_a = _1] >> *catchall >
          raw_token(T_END) > tok.ident[_pass = (_1 == _a)] ;

      via = raw_token(T_VIA) > tok.ident[_a = _1] >>
         *(catchall | raw_token(T_LAYER) | raw_token(T_FOREIGN)) >>
         raw_token(T_END) > tok.ident[_pass = (_1 == _a)] ;

      viarule = raw_token(T_VIARULE) > tok.ident[_a = _1] >>
	*(catchall | raw_token(T_LAYER) | raw_token(T_BY) | raw_token(T_VIA)) >>
        raw_token(T_END) > tok.ident[_pass = (_1 == _a)] ;

      // BOZO review case sensitivity of these keywords and correct if necessary
      // BOZO this can appear in DEF as well; move to common LEF/DEF header
      // but how? this is internal to class.  Oh hey, inheritance...?

      // one-symbol site classes
      siteclass = (raw_token(T_COVER) >> attr(SITECLASS_COVER)) |
	          (raw_token(T_RING) >> attr(SITECLASS_RING)) |
	          (raw_token(T_PAD) >> attr(SITECLASS_PAD)) |
	          (raw_token(T_BLOCK) >> attr(SITECLASS_BLOCK)) |
	          (raw_token(T_CORE) >> attr(SITECLASS_CORE));
      // classes distinguished by second symbol (all ENDCAP and some PAD, CORE)
      padclass =  (raw_token(T_INPUT) >> attr(SITECLASS_PAD_INPUT)) |
	          (raw_token(T_OUTPUT) >> attr(SITECLASS_PAD_OUTPUT)) |
	          (raw_token(T_SPACER) >> attr(SITECLASS_PAD_SPACER)) |
	          (raw_token(T_INOUT) >> attr(SITECLASS_PAD_INOUT)) |
	          (raw_token(T_POWER) >> attr(SITECLASS_PAD_POWER));
      coreclass = (raw_token(T_FEEDTHRU) >> attr(SITECLASS_CORE_FEEDTHRU)) |
	          (raw_token(T_TIEHIGH) >> attr(SITECLASS_CORE_TIEHIGH)) |
	          (raw_token(T_TIELOW) >> attr(SITECLASS_CORE_TIELOW));
      endcapclass = (raw_token(T_PRE) >> attr(SITECLASS_ENDCAP_PRE)) ||
	          (raw_token(T_POST) >> attr(SITECLASS_ENDCAP_POST)) ||
	          (raw_token(T_TOPLEFT) >> attr(SITECLASS_ENDCAP_TOPLEFT)) ||
	          (raw_token(T_TOPRIGHT) >> attr(SITECLASS_ENDCAP_TOPRIGHT)) ||
	          (raw_token(T_BOTTOMLEFT) >> attr(SITECLASS_ENDCAP_BOTTOMLEFT)) ||
	          (raw_token(T_BOTTOMRIGHT) >> attr(SITECLASS_ENDCAP_BOTTOMRIGHT));

      classrule = raw_token(T_CLASS) >
	          ((raw_token(T_PAD) >> padclass) |    // three kinds of two-keyword sites
		   (raw_token(T_CORE) >> coreclass) |
		   (raw_token(T_ENDCAP) > endcapclass) |
		   siteclass) > ';' ;               // try single keyword sites last

      sitesym = (raw_token(T_X) >> attr(SITESYM_X)) |
	        (raw_token(T_Y) >> attr(SITESYM_Y)) |
	        (raw_token(T_R90) >> attr(SITESYM_R90));

      // BOZO add encountered sites to symbol table for checking
      // BOZO need a better way to handle case-insensitive keywords
      siterule %= raw_token(T_SITE) > tok.ident[_a = _1] >
	classrule >
	-(raw_token(T_SYMMETRY) > +sitesym > ';' ) >
	// BOZO use lefextent here
	raw_token(T_SIZE) > tok.double_ > raw_token(T_BY) > tok.double_ > ';' >
	// BOZO optional ROWPATTERN goes here
        raw_token(T_END) > omit[tok.ident[_pass = (_1 == _a)]] ;

      spacing = raw_token(T_SPACING) >> *catchall >> raw_token(T_END) > raw_token(T_SPACING) ;
      units = raw_token(T_UNITS) >> *catchall >> raw_token(T_END) >> raw_token(T_UNITS) ;

      // macro properties
      point = tok.double_ >> tok.double_ ;
      foreign = raw_token(T_FOREIGN) > tok.ident > point > ';' ;
      origin  = raw_token(T_ORIGIN) > point > ';' ;
      macrosymmetry = raw_token(T_SYMMETRY) > *sitesym > ';' ;
      macrosize = raw_token(T_SIZE) > tok.double_ > no_case[raw_token(T_BY)] > tok.double_ > ';' ;

      // macro elements
      pin_geom = raw_token(T_PORT) >> *(raw_token(T_LAYER) > tok.ident > ';' > +(tok.ident | tok.double_ | ';')) > raw_token(T_END) ;
      pin = raw_token(T_PIN) > tok.ident[_a = _1] >>
	*(catchall | padclass | pin_geom) >>
        raw_token(T_END) >> tok.ident[_pass = (_1 == _a)] ;
      obs = raw_token(T_OBS) > *(catchall | raw_token(T_LAYER)) > raw_token(T_END) ;

      macro = raw_token(T_MACRO) > tok.ident[_a = _1, at_c<0>(_val) = _1] >
	// these can be in any order, but the caret checks that at least one is present
	// SIZE is required, so I guess that makes it work
	// have to use semantic actions here because the () grouped stuff becomes its own sequence
	(classrule[at_c<1>(_val) = _1] ^
	 foreign[at_c<2>(_val) = _1] ^
	 origin[at_c<3>(_val) = _1] ^
	 macrosize[at_c<4>(_val) = _1] ^
	 macrosymmetry[at_c<5>(_val) = _1] ^
	 (raw_token(T_SITE) > tok.ident[at_c<6>(_val) = _1] > ';')) >
	// BOZO not stored for now:
	*pin >
	-obs >
        raw_token(T_END) > tok.ident[_pass = (_1 == _a)] ;

      layerdef = raw_token(T_LAYER) > tok.ident[_a = _1, at_c<0>(_val) = _1]
        > raw_token(T_TYPE)
        >   (raw_token(T_CUT)        [at_c<1>(_val) = Cut]
           | raw_token(T_MASTERSLICE)[at_c<1>(_val) = Masterslice]
           | raw_token(T_ROUTING)    [at_c<1>(_val) = Routing]
           | raw_token(T_IMPLANT)    [at_c<1>(_val) = Implant]
           | raw_token(T_OVERLAP)    [at_c<1>(_val) = Overlap])
        > *catchall
        > raw_token(T_END) > tok.ident[_pass = (_a == _1)] ;

      // define some major elements

      // strangely a LEF file does not begin with anything distinctive... just goes right into the data
      // BOZO siterule should store site name into symtab
      // evidently EVERY statement is optional
      lef_file = *(casesens | units | spacing | via | viarule |
                   version_stmt[at_c<3>(_val) = _1] |
                   layerdef [push_back(at_c<2>(_val), _1)] |
		   siterule [push_back(at_c<0>(_val), _1)] |
		   macro [push_back(at_c<1>(_val), _1)])
	>> -(raw_token(T_END) > raw_token(T_LIBRARY)) ;

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

  boost::spirit::qi::rule<Iterator, double()> version_stmt;

  boost::spirit::qi::rule<Iterator, leflayer(),
                          boost::spirit::qi::locals<std::string> > layerdef;

  // for now macro synthesizes a string AND has a local variable (for checking END statement)
  boost::spirit::qi::rule<Iterator, lefmacro(),
			  boost::spirit::qi::locals<std::string> > macro;

  // Some keywords turn directly into enums in the database
  boost::spirit::qi::rule<Iterator, SiteClass()> siteclass, padclass, coreclass, endcapclass;
  boost::spirit::qi::rule<Iterator, SiteSymmetry()> sitesym;
  boost::spirit::qi::rule<Iterator, SiteClass()> classrule;
  // site rule uses local var (to check END statement) and synthesizes a Site
  boost::spirit::qi::rule<Iterator, Site(),
			  boost::spirit::qi::locals<std::string> > siterule;
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
				  boost::spirit::qi::locals<std::string> > Unparsed;
  Unparsed layer, via, viarule, pin, pin_geom, catchall;

  // stuff inside macros I don't handle
  boost::spirit::qi::rule<Iterator > unparsed_macro_stuff;

  // The LEF file as a whole
  boost::spirit::qi::rule<Iterator, lef()> lef_file;

  // error handler
  boost::phoenix::function<error_info_impl> error_info;
};

}

#endif
