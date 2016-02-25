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

// This file contains result types for LEF parsing

#if !defined(EDASKEL_LEF_TYPES)
#define EDASKEL_LEF_TYPES

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>
#include <string>

#include "../db/globals.h"
#include "deftypes.h"

using namespace EDASkel;

namespace LefParse {

struct leflayer {
  std::string name;
  LayerTypes type;
};

// I'm naming this differently from the other types in here. Is that good?
// should this be global?
struct Site {
  std::string name;
  SiteClass class_;
  boost::optional<std::vector<SiteSymmetry> > symmetry;
  double width, height;
};

struct lefpoint {
  double x, y;
};

struct lefextent {
  double width, height;
};

struct lefforeign {
  std::string name;
  lefpoint pt;
};

struct lefmacro {
  std::string name;
  boost::optional<SiteClass> class_;
  boost::optional<lefforeign> foreign;
  boost::optional<lefpoint> origin;
  boost::optional<lefextent> size;
  boost::optional<std::vector<SiteSymmetry> > symmetry;
  boost::optional<std::string> site;
};  

struct lef {
  double version;
  std::vector<Site> sites;
  std::vector<lefmacro> macros;
  std::vector<leflayer> layers;
};

}

// BOOST_FUSION_ADAPT_STRUCT needs to be outside the namespace for some reason

BOOST_FUSION_ADAPT_STRUCT(
  LefParse::leflayer,
  (std::string, name)
  (LayerTypes, type)
)

BOOST_FUSION_ADAPT_STRUCT(
  LefParse::Site,
  (std::string, name)
  (SiteClass, class_)
  (boost::optional<std::vector<SiteSymmetry> >, symmetry)
  (double, width)
  (double, height)
)

BOOST_FUSION_ADAPT_STRUCT(
  LefParse::lefpoint,
  (double, x)
  (double, y)
)

BOOST_FUSION_ADAPT_STRUCT(
  LefParse::lefextent,
  (double, width)
  (double, height)
)

BOOST_FUSION_ADAPT_STRUCT(
  LefParse::lefforeign,
  (std::string, name)
  (LefParse::lefpoint, pt)
)

BOOST_FUSION_ADAPT_STRUCT(
  LefParse::lefmacro,
  (std::string, name)
  (boost::optional<SiteClass>, class_)
  (boost::optional<LefParse::lefforeign>, foreign)
  (boost::optional<LefParse::lefpoint>, origin)
  (boost::optional<LefParse::lefextent>, size)
  (boost::optional<std::vector<SiteSymmetry> >, symmetry)
  (boost::optional<std::string>, site)
)

BOOST_FUSION_ADAPT_STRUCT(
  LefParse::lef,
  (std::vector<LefParse::Site>, sites)
  (std::vector<LefParse::lefmacro>, macros)
  (std::vector<LefParse::leflayer>, layers)
  (double, version)
)

#endif
