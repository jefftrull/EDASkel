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

// This file contains result types for DEF parsing
#if !defined(EDASKEL_DEF_TYPES)
#define EDASKEL_DEF_TYPES


#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/optional.hpp>
#include <string>

namespace DefParse {

struct defpoint {
  int x, y;
};

struct defrect {
  defpoint ll, ur;
};

struct defplcinfo {
  std::string plcfix;       // BOZO make boolean "fixed" or something
  defpoint origin;
  std::string orient;       // BOZO make this enum?
};

struct defcomponent {
  std::string name;
  std::string celltype;
  boost::optional<defplcinfo> placement;
};

// BOOST_FUSION_ADAPT_STRUCT is a macro and will be confused by this embedded comma unless:
typedef std::pair<int, int> IntPair;
struct siterepeat {
  int xrepeat;
  int yrepeat;
  boost::optional<IntPair> step;
};

struct rowsite_b {
  std::string sitename;
  int origx;
  int origy;
  std::string orient;
  boost::optional<siterepeat> repeat;
};

struct rowsite {
  boost::optional<std::string> rowname;
  rowsite_b body;
};

struct def {
  std::string name;
  double version;
  defrect diearea;
  int dbupermicron;
  std::vector<defcomponent> components;
  std::vector<rowsite> rows;
};

}

BOOST_FUSION_ADAPT_STRUCT(
  DefParse::defpoint,
  (int, x)
  (int, y)
)

BOOST_FUSION_ADAPT_STRUCT(
  DefParse::defrect,
  (DefParse::defpoint, ll)
  (DefParse::defpoint, ur)
)

BOOST_FUSION_ADAPT_STRUCT(
  DefParse::defplcinfo,
  (std::string, plcfix)
  (DefParse::defpoint, origin)
  (std::string, orient)
)

BOOST_FUSION_ADAPT_STRUCT(
  DefParse::defcomponent,
  (std::string, name)
  (std::string, celltype)
  (boost::optional<DefParse::defplcinfo>, placement)
)

BOOST_FUSION_ADAPT_STRUCT(
  DefParse::siterepeat,
  (int, xrepeat)
  (int, yrepeat)
  (boost::optional<DefParse::IntPair>, step)
)

BOOST_FUSION_ADAPT_STRUCT(
  DefParse::rowsite_b,
  (std::string, sitename)
  (int, origx)
  (int, origy)
  (std::string, orient)
  (boost::optional<DefParse::siterepeat>, repeat)
)

BOOST_FUSION_ADAPT_STRUCT(
  DefParse::rowsite,
  (boost::optional<std::string>, rowname)
  (DefParse::rowsite_b, body)
)

BOOST_FUSION_ADAPT_STRUCT(
  DefParse::def,
  (std::string, name)
  (double, version)
  (DefParse::defrect, diearea)
  (int, dbupermicron)
  (std::vector<DefParse::defcomponent>, components)
  (std::vector<DefParse::rowsite>, rows)
)

#endif
