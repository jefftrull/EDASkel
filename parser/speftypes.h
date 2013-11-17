// Types (AST, basically) used in parsing SPEF
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

#ifndef PARSER_SPEFTYPES_H
#define PARSER_SPEFTYPES_H

#define FUSION_MAX_VECTOR_SIZE 15

#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/resistance.hpp>
#include <boost/units/systems/si/capacitance.hpp>
#include <boost/units/systems/si/inductance.hpp>
#include <boost/units/systems/si/io.hpp>

#include <string>
#include <unordered_map>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

namespace EDASkel {
  namespace SpefParse {
    using namespace boost::units;
    typedef quantity<si::time, double>         time_units_t;
    typedef quantity<si::resistance, double>   resistance_units_t;
    typedef quantity<si::capacitance, double>  capacitance_units_t;
    typedef quantity<si::inductance, double>   inductance_units_t;
    typedef std::map<std::string, std::string> design_flow_map_t;
    typedef unsigned name_map_index_t;
    typedef std::unordered_map<name_map_index_t, std::string> name_map_t;
    struct spef {
      std::string standard;          // i.e. IEEE 1481-1998 - the only one we support for now
      std::string name;
      boost::posix_time::ptime date; // a date *and* time, actually
      std::string         vendor, program, version;
      design_flow_map_t   design_flow;
      time_units_t        t_unit;     // multiplier for all duration values
      capacitance_units_t c_unit;     // capacitance
      resistance_units_t  r_unit;     // and for resistance
      inductance_units_t  l_unit;     // and inductance
      name_map_t          name_map;
    };
  }
}

BOOST_FUSION_ADAPT_STRUCT(
  EDASkel::SpefParse::spef,
  (std::string,                             standard)
  (std::string,                             name)
  (std::string,                             vendor)
  (std::string,                             program)
  (std::string,                             version)
  (EDASkel::SpefParse::design_flow_map_t,   design_flow)
  (EDASkel::SpefParse::time_units_t,        t_unit)
  (EDASkel::SpefParse::capacitance_units_t, c_unit)
  (EDASkel::SpefParse::resistance_units_t,  r_unit)
  (EDASkel::SpefParse::inductance_units_t,  l_unit)
  (EDASkel::SpefParse::name_map_t,          name_map)
)

#endif // PARSER_SPEFTYPES_H
