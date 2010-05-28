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

// This file contains result types for DEF parsing
#if !defined(EDASKEL_DEF_TYPES)
#define EDASKEL_DEF_TYPES


#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>
#include <string>

struct defpoint {
  int x, y;
};

BOOST_FUSION_ADAPT_STRUCT(
  defpoint,
  (int, x)
  (int, y)
)

struct defrect {
  defpoint ll, ur;
};

BOOST_FUSION_ADAPT_STRUCT(
  defrect,
  (defpoint, ll)
  (defpoint, ur)
)

struct defplcinfo {
  std::string plcfix;       // BOZO make boolean "fixed" or something
  defpoint origin;
  std::string orient;       // BOZO make this enum
};

BOOST_FUSION_ADAPT_STRUCT(
  defplcinfo,
  (std::string, plcfix)
  (defpoint, origin)
  (std::string, orient)
)

struct defcomponent {
  std::string name;
  std::string celltype;
  boost::optional<defplcinfo> placement;
};

BOOST_FUSION_ADAPT_STRUCT(
  defcomponent,
  (std::string, name)
  (std::string, celltype)
  (boost::optional<defplcinfo>, placement)
)

struct def {
  std::string name;
  double version;
  defrect diearea;
  std::vector<defcomponent> components;
};

BOOST_FUSION_ADAPT_STRUCT(
  def,
  (std::string, name)
  (double, version)
  (defrect, diearea)
  (std::vector<defcomponent>, components)
)

#endif
