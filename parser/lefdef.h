// A Boost Spirit-based LEF/DEF parser, part of EDASkel, a sample EDA app
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

// This file contains common definitions for the LEF and DEF parsers

#if !defined(EDASKEL_LEFDEF)
#define EDASKEL_LEFDEF

#include <boost/spirit/include/lex_lexertl.hpp>
#include <boost/spirit/include/lex_lexertl_position_token.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

namespace EDASkel {

// A skip parser for LEF/DEF comments and spaces
// adapted from a presentation at Boostcon 2010 by Michael Caisse
template <typename Iterator>
struct lefdefskipper : boost::spirit::qi::grammar< Iterator >
{
 lefdefskipper() : lefdefskipper::base_type(skip_it)
    {
      using namespace boost::spirit::qi;

      comment = '#' >> *( char_ - eol ) >> eol ;
      skip_it = comment | space ;
    }
  boost::spirit::qi::rule<Iterator> skip_it;
  boost::spirit::qi::rule<Iterator> comment;
};

  // typedef for stream iterator we will use
  typedef boost::spirit::istream_iterator LefDefIter;
  // lexer needs the iterator type and a list of token attribute types
  typedef boost::spirit::lex::lexertl::position_token<LefDefIter,
    boost::mpl::vector<int, double, std::string> > LefDefToken;
  typedef boost::spirit::lex::lexertl::actor_lexer<LefDefToken> LefDefLexer;

}

#endif
