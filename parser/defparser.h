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

#include <boost/spirit/include/lex_lexertl.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/lexical_cast.hpp>

#include "deftypes.h"

namespace DefParse {

// DEF tokens
enum TokenIds {
  T_VERSION = 1000,        // can't start at 0 == EOF
  T_DIEAREA, T_WEIGHT, T_SOURCE, T_DIST, T_NETLIST,
  T_USER, T_TIMING, T_COMPONENTS, T_END, T_DO, T_BY,
  T_STEP, T_ROW, T_SITE, T_UNITS, T_DISTANCE, T_MICRONS,
  T_TRACKS, T_GCELLGRID, T_DESIGN, T_PINS, T_NETS,
  T_SPECIALNETS, T_VIAS, T_PLACED, T_FIXED, T_ANY
};


template <typename Lexer>
struct DefTokens : boost::spirit::lex::lexer<Lexer>
{
  DefTokens()
    : history_("^HISTORY[^;]*;"),
      double_(R"(-?[0-9]+\.[0-9]+)"),
      int_("-?[0-9]+")
  {
    namespace lex = boost::spirit::lex;
    // for lex semantic actions
    using lex::_val;
    using lex::_start;
    using lex::_end;
    using namespace boost::phoenix::local_names;
    using boost::phoenix::construct;
    using boost::phoenix::let;
    using boost::phoenix::begin;
    using boost::phoenix::end;
    namespace phx = boost::phoenix;

    // identifiers
    // we cannot distinguish in the lexer between the different valid instance, cell, and design names
    // because the lexer has no context and considers each pattern in order
    // so we must do the most general possibility and, perhaps, check validity in the parser

    // some identifiers (esp instance names) may be as complex as:
    // letter followed by letters, numbers, underscore, hyphen, square brackets, slashes (hierarchy)
    // with potentially embedded, quoted bracketed numbers, and optionally a final unquoted bracketed number
    nonkwd_ = R"([a-zA-Z]([a-zA-Z0-9_/]|-|(\\\[[0-9+]\\\]))*(\[[0-9]+\])?)";
    
    this->self +=
      // history is highest priority - may contain keywords!
      // cleaner, but gives syntax error:
//      history_ [ _val = construct<std::string>(_val, 8, size(_val)-9) ]
        history_ [ let(_a = construct<std::string>(_start, _end))
                  // have to create a string to do this because _start/_end are istream iterators
                  // and cannot be indexed (random access)
                  [ _val = construct<std::string>(begin(_a)+8, end(_a) -1) ]
          ]
      | lex::string("VERSION", T_VERSION)
      | lex::string("DIEAREA", T_DIEAREA)
      | lex::string("WEIGHT", T_WEIGHT)
      | lex::string("SOURCE", T_SOURCE)
      | lex::string("DIST", T_DIST)
      | lex::string("NETLIST", T_NETLIST)
      | lex::string("USER", T_USER)
      | lex::string("TIMING", T_TIMING)
      | lex::string("COMPONENTS", T_COMPONENTS)
      | lex::string("END", T_END)
      | lex::string("DO", T_DO)
      | lex::string("BY", T_BY)
      | lex::string("STEP", T_STEP)
      | lex::string("ROW", T_ROW)
      | lex::string("SITE", T_SITE)
      | lex::string("UNITS", T_UNITS)
      | lex::string("DISTANCE", T_DISTANCE)
      | lex::string("MICRONS", T_MICRONS)
      | lex::string("TRACKS", T_TRACKS)
      | lex::string("GCELLGRID", T_GCELLGRID)
      | lex::string("DESIGN", T_DESIGN)
      | lex::string("PINS", T_PINS)
      | lex::string("NETS", T_NETS)
      | lex::string("SPECIALNETS", T_SPECIALNETS)
      | lex::string("VIAS", T_VIAS)
      | lex::string("PLACED", T_PLACED)
      | lex::string("FIXED", T_FIXED)
      | nonkwd_
      | double_ | int_
      | '+' | '-' | '(' | ')' | ';'
      ;

    // whitespace and comments
    this->self += lex::string(R"(([ \t\n]+)|(#[^\n]*\n))")
        [
            lex::_pass = lex::pass_flags::pass_ignore
        ];

    // catchall (mostly for error processing)
    this->self += lex::string(".", T_ANY);
  }

  // token definitions
  // attribute-less tokens are all done above by tokenid per hkaiser recommendation

  // string attribute tokens (different kinds of identifiers, mostly)
  boost::spirit::lex::token_def<std::string> nonkwd_, history_;
  // numbers
  boost::spirit::lex::token_def<double> double_;
  boost::spirit::lex::token_def<int> int_;

};

// a starter DEF grammar
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

// first define some major subcomponents in their own "grammar"

// Components
typedef boost::spirit::qi::symbols<char, std::string> comp_symtab_t;  // for reference checking

template <typename Iterator, typename Lexer>
struct comp_parser : boost::spirit::qi::grammar<Iterator,
                                                boost::spirit::locals<std::string>,
                                                defcomponent()>
{
  template <typename TokenDef>
    comp_parser(TokenDef const& tok, comp_symtab_t& compsym) :
   comp_parser::base_type(comp), comp_symtab(compsym)
  {
    using namespace boost::spirit::qi;
    using boost::phoenix::bind;
    using boost::phoenix::at_c;
    using boost::spirit::_1;
    using boost::spirit::_val;

    point = '(' >> tok.int_ >> tok.int_ >> ')' ;       // points are parenthesized pairs, no comma

    orient %= tok.nonkwd_[_pass = ((_1 == "N")  || (_1 == "S") ||  (_1 == "E")  || (_1 == "W") ||
                                   (_1 == "FN") || (_1 == "FS") || (_1 == "FE") || (_1 == "FW"))] ;

    plcinfo = '+' >> ((token(T_FIXED) > point > orient) |
                       (token(T_PLACED) > point > orient)) ;    // location and orientation

    weight = '+' >> (raw_token(T_WEIGHT) > tok.int_) ;

    source = '+' >> (raw_token(T_SOURCE) > (raw_token(T_DIST) | raw_token(T_NETLIST) | raw_token(T_USER) | raw_token(T_TIMING))) ;

    // components required instance name and celltype; optional any (or none) of placement and weight, in any order:
    comp %=  '-' > tok.nonkwd_[_a = _1] > tok.nonkwd_[bind(comp_symtab.add, _a, _1)]
           > (plcinfo ^ omit[weight] ^ omit[source] ^ eps) > ';' ;

    BOOST_SPIRIT_DEBUG_NODE(point);
    BOOST_SPIRIT_DEBUG_NODE(weight);
    BOOST_SPIRIT_DEBUG_NODE(source);
    BOOST_SPIRIT_DEBUG_NODE(comp);

    on_error<fail>
    (
      comp
    , std::cerr << error_info(boost::spirit::_3, boost::spirit::_4) << std::endl
    );

  }

  // top level COMPONENT
  boost::spirit::qi::rule<Iterator, boost::spirit::locals<std::string>, defcomponent()> comp;

  // Symbol table storage for components
  comp_symtab_t& comp_symtab;

  // points "( x y )" produces defpoint structs (see deftypes.h)
  boost::spirit::qi::rule<Iterator, defpoint()> point;

  // orientation synthesizes a string
  boost::spirit::qi::rule<Iterator, std::string()> orient;

  // optional placement info (placed vs. fixed, location, orientation)
  boost::spirit::qi::rule<Iterator, defplcinfo()> plcinfo;

  // WEIGHT/SOURCE - presently no attributes or locals (parsed but not stored)
  boost::spirit::qi::rule<Iterator> weight, source;

  // error handler
  boost::phoenix::function<error_info_impl> error_info;
};   

// A NET statement
template <typename Iterator, typename Lexer>
struct net_parser : boost::spirit::qi::grammar<Iterator, defnet()>
{
  template <typename TokenDef>
  net_parser(TokenDef const& tok, comp_symtab_t const& compsym) :
  net_parser::base_type(net), comp_symtab(compsym)
  {
    using namespace boost::spirit::qi;

    // symbol tables cannot be used directly in parsers taking lexer tokens
    // however, we can get the token's string value and then do a lookup
    // we will do this for the component names supplied in connections

    // disambiguate call to symbols::find
    typedef typename comp_symtab_t::value_type const* (comp_symtab_t::*findfn_t)(std::string const&) const;

    using boost::spirit::_1;
    using boost::spirit::_a;
    namespace phx = boost::phoenix;
    connection %= '(' > tok.nonkwd_[_a = _1] // component name
                > eps[_pass = phx::bind(static_cast<findfn_t>(&comp_symtab_t::find),
                                        phx::cref(comp_symtab), _a)]
                > tok.nonkwd_ > ')' ;

    net = '-' > tok.nonkwd_ > *connection > ';' ;
  }

  boost::spirit::qi::rule<Iterator,
                          boost::spirit::locals<std::string>,
                          defconnection()> connection;
  boost::spirit::qi::rule<Iterator, defnet()> net;

  comp_symtab_t const& comp_symtab;
};   


template <typename Iterator, typename Lexer>
struct defparser : boost::spirit::qi::grammar<Iterator, def()>
{

  template <typename TokenDef>
    defparser(TokenDef const& tok) : defparser::base_type(def_file),
                                     component(tok, compsym), net(tok, compsym)
    {
      using namespace boost::spirit::qi;
      using boost::spirit::_1;                  // access attributes for component count check
      using boost::phoenix::push_back;          // to store results in containers
      using boost::phoenix::at_c;               // to refer to pieces of wrapped structs

      // top-level elements in a DEF file
      version_stmt = raw_token(T_VERSION) > tok.double_ > ';' ;

      point = '(' >> tok.int_ >> tok.int_ >> ')' ;       // points are parenthesized pairs, no comma
      rect = point >> point ;                            // rects are just two points in a row
      diearea_stmt = raw_token(T_DIEAREA) > rect > ';' ;

      // define some major elements

      // define repeat rules for each element section
      // I tried to make this generic but failed.  You can pass in the rule as an inherited attribute,
      // but the type of the parent rule needs to be different for each since the synthesized attribute is
      // different (i.e., std::vector<attribute_of_child_rule>)
      comps_section %= raw_token(T_COMPONENTS) > omit[tok.int_[_a = _1]] > ';' >  // remember count in a local variable
	               repeat(_a)[component] >                           // expect that many copies of the supplied rule
                       raw_token(T_END) > raw_token(T_COMPONENTS) ;      // END followed by the section name again

      // Similarly, other counted (but currently unparsed) stuff:
      vias_section %= raw_token(T_VIAS) > omit[tok.int_[_a = _1]] > ';' > 
	              repeat(_a)['-' > *~char_(';') > ';'] >
                      raw_token(T_END) > raw_token(T_VIAS) ;
      nets_section %= raw_token(T_NETS) > omit[tok.int_[_a = _1]] > ';' > 
	              repeat(_a)[net] >
                      raw_token(T_END) > raw_token(T_NETS) ;
      specialnets_section %= raw_token(T_SPECIALNETS) > omit[tok.int_[_a = _1]] > ';' > 
	                     repeat(_a)['-' > *~char_(';') > ';'] >
                             raw_token(T_END) > raw_token(T_SPECIALNETS) ;
      pins_section %= raw_token(T_PINS) > omit[tok.int_[_a = _1]] > ';' > 
	              repeat(_a)['-' > *~char_(';') > ';'] >
                      raw_token(T_END) > raw_token(T_PINS) ;

      // My copy of the LEF/DEF reference does not show this SITE command as valid for DEF yet my example data does...
      // The example data's syntax is very similar to that defined for ROW, so I'll combine them
      siterpt_stmt = raw_token(T_DO) > tok.int_ > raw_token(T_BY) > tok.int_ > -(raw_token(T_STEP) > tok.int_ > tok.int_) ;
      orient %= tok.nonkwd_[_pass = ((_1 == "N")  || (_1 == "S") ||  (_1 == "E")  || (_1 == "W") ||
                                     (_1 == "FN") || (_1 == "FS") || (_1 == "FE") || (_1 == "FW"))] ;
      rowsite_stmt = ((raw_token(T_ROW) > tok.nonkwd_) | raw_token(T_SITE)) > tok.nonkwd_ > tok.int_ > tok.int_ > orient >
	             -siterpt_stmt > ';' ;

      dbu = raw_token(T_UNITS) > raw_token(T_DISTANCE) > raw_token(T_MICRONS) > tok.int_ > ';' ;

      // This parser only handles components and a couple of misc. statements
      // here's a catchall parser to discard all other data
      tracks_stmt = raw_token(T_TRACKS) >> *(tok.nonkwd_ | raw_token(T_DO) | raw_token(T_STEP) | tok.int_) > ';' ;
      gcellgrid_stmt = raw_token(T_GCELLGRID) > *(tok.nonkwd_ | raw_token(T_DO) | raw_token(T_STEP) | tok.int_) > ';' ;

      unparsed = vias_section | specialnets_section | pins_section | tracks_stmt | gcellgrid_stmt | history_stmt ;

      def_file = raw_token(T_DESIGN) > tok.nonkwd_[at_c<0>(_val) = _1] > ';' >
                 *(version_stmt[at_c<1>(_val) = _1] |
		   diearea_stmt[at_c<2>(_val) = _1] |
		   dbu[at_c<3>(_val) = _1] |
      	           comps_section[at_c<4>(_val) = _1] |
      	           nets_section[at_c<5>(_val) = _1] |
		   rowsite_stmt[push_back(at_c<6>(_val), _1)] |
                   tok.history_[push_back(at_c<7>(_val), _1)] |
		   unparsed) >
	raw_token(T_END) > raw_token(T_DESIGN) ;

      // Debugging assistance

      BOOST_SPIRIT_DEBUG_NODE(version_stmt);
      BOOST_SPIRIT_DEBUG_NODE(diearea_stmt);
      BOOST_SPIRIT_DEBUG_NODE(comps_section);
      BOOST_SPIRIT_DEBUG_NODE(tracks_stmt);

      on_error<fail>
      (
        def_file
      , std::cerr << error_info(boost::spirit::_3, boost::spirit::_4) << std::endl
      );
    }

  // a rule representing the entire COMPONENTS section
  boost::spirit::qi::rule<Iterator, std::vector<defcomponent>(), boost::spirit::qi::locals<int> > comps_section;

  // a symbol table holding all the components that were found
  comp_symtab_t compsym;

  // a single instance within the COMPONENTS section (name, celltype, placement)
  comp_parser<Iterator, Lexer> component;

  // The entire NETS section
  boost::spirit::qi::rule<Iterator, std::vector<defnet>(), boost::spirit::qi::locals<int> > nets_section;

  // a single net within the NETS section
  net_parser<Iterator, Lexer> net;

  // VERSION takes no parameters (a.k.a. "inherited attributes") and synthesizes a double for its attribute
  boost::spirit::qi::rule<Iterator, double()> version_stmt;

  // points "( x y )" produces defpoint structs (see deftypes.h)
  boost::spirit::qi::rule<Iterator, defpoint()> point;
  // rects "( llx lly ) ( urx ury )" synthesize defrect structs containing two points
  boost::spirit::qi::rule<Iterator, defrect()> rect, diearea_stmt;

  typedef boost::spirit::qi::rule<Iterator, std::string()> StringRule;

  // Site (or named row) statement
  boost::spirit::qi::rule<Iterator, std::string()> orient;
  boost::spirit::qi::rule<Iterator, siterepeat()> siterpt_stmt;
  boost::spirit::qi::rule<Iterator, rowsite()> rowsite_stmt;

  boost::spirit::qi::rule<Iterator, int()> dbu;

  // a catchall rule for everything I don't (yet) parse.  No attribute synthesized.
  boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<std::string> > tracks_stmt, gcellgrid_stmt, history_stmt;
  // stuff I count, but don't parse
  boost::spirit::qi::rule<Iterator, boost::spirit::qi::locals<int> > vias_section, specialnets_section, pins_section;
  // top level unparsed section - no attribute, no local variables
  boost::spirit::qi::rule<Iterator> unparsed;

  // The DEF file as a whole
  boost::spirit::qi::rule<Iterator, def()> def_file;

  // error handler
  boost::phoenix::function<error_info_impl> error_info;
};

}

#endif
