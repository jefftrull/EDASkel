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

// Compiling the LEF and DEF parsers together is very slow on my system
// (large memory requirement), and many files (esp. unit tests) depend on
// them, so this is a separate "compilation unit" that uses the parsers.
// I will make other files depend on this one at link time.

#include <boost/spirit/include/support_multi_pass.hpp>

#include "lefparser.h"
#include "defparser.h"
using namespace LefParse;
using namespace DefParse;
using namespace EDASkel;

namespace EDASkel {
  typedef std::string::const_iterator LefDefStringIter;
  defparser<LefDefStringIter> defStringParser;
  lefparser<LefDefStringIter> lefStringParser;
  lefdefskipper<LefDefStringIter> lefdefStringSkipper;

  // use a special adapter for files (stream iterators are "Input", not "Forward"
  // iterators, and so need some extra help to handle lookahead)
  typedef boost::spirit::multi_pass<std::istreambuf_iterator<char> > LefDefFileIter;
  defparser<LefDefFileIter> defFileParser;
  lefparser<LefDefFileIter> lefFileParser;
  lefdefskipper<LefDefFileIter> lefdefFileSkipper;
}

