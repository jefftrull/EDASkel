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

// This file contains result types for LEF parsing

#if !defined(EDASKEL_LEF_TYPES)
#define EDASKEL_LEF_TYPES

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>
#include <string>

#include "deftypes.h"

// these ought to move into some LEF/DEF common file
enum SiteSymmetry { SITESYM_X, SITESYM_Y, SITESYM_R90 };
enum SiteClass { SITECLASS_COVER, SITECLASS_RING, SITECLASS_BLOCK,
		 SITECLASS_PAD, SITECLASS_PAD_INPUT, SITECLASS_PAD_OUTPUT, SITECLASS_PAD_INOUT,
		 SITECLASS_PAD_POWER, SITECLASS_PAD_SPACER,
		 SITECLASS_CORE, SITECLASS_CORE_FEEDTHRU, SITECLASS_CORE_TIEHIGH, SITECLASS_CORE_TIELOW,
		 SITECLASS_ENDCAP_PRE, SITECLASS_ENDCAP_POST, SITECLASS_ENDCAP_TOPLEFT,
		 SITECLASS_ENDCAP_TOPRIGHT, SITECLASS_ENDCAP_BOTTOMLEFT, SITECLASS_ENDCAP_BOTTOMRIGHT };
struct Site {
  std::string name;
  SiteClass class_;
  boost::optional<SiteSymmetry> symmetry;
  float width, height;
};

BOOST_FUSION_ADAPT_STRUCT(
  Site,
  (std::string, name)
  (SiteClass, class_)
  (boost::optional<SiteSymmetry>, symmetry)
  (float, width)
  (float, height)
)

struct lefpoint {
  float x, y;
};

BOOST_FUSION_ADAPT_STRUCT(
  lefpoint,
  (float, x)
  (float, y)
)

struct lefextent {
  float width, height;
};

BOOST_FUSION_ADAPT_STRUCT(
  lefextent,
  (float, width)
  (float, height)
)

struct lefforeign {
  std::string name;
  lefpoint pt;
};

BOOST_FUSION_ADAPT_STRUCT(
  lefforeign,
  (std::string, name)
  (lefpoint, pt)
)

struct lefmacro {
  std::string name;
  boost::optional<SiteClass> class_;
  boost::optional<lefforeign> foreign;
  boost::optional<lefpoint> origin;
  boost::optional<lefextent> size;
  boost::optional<std::vector<SiteSymmetry> > symmetry;
  boost::optional<std::string> site;
};  

BOOST_FUSION_ADAPT_STRUCT(
  lefmacro,
  (std::string, name)
  (boost::optional<SiteClass>, class_)
  (boost::optional<lefforeign>, foreign)
  (boost::optional<lefpoint>, origin)
  (boost::optional<lefextent>, size)
  (boost::optional<std::vector<SiteSymmetry> >, symmetry)
  (boost::optional<std::string>, site)
)

struct lef {
  std::vector<Site> sites;
  std::vector<lefmacro> macros;
};

BOOST_FUSION_ADAPT_STRUCT(
  lef,
  (std::vector<Site>, sites)
  (std::vector<lefmacro>, macros)
)

#endif
