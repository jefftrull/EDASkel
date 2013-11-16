#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <boost/units/systems/si/time.hpp>
#include <boost/units/systems/si/resistance.hpp>
#include <boost/units/systems/si/capacitance.hpp>
#include <boost/units/systems/si/inductance.hpp>
#include <boost/units/systems/si/io.hpp>

#include <string>

#include <boost/fusion/include/adapt_struct.hpp>

#ifndef PARSER_SPEFTYPES_H
#define PARSER_SPEFTYPES_H

namespace EDASkel {
  namespace SpefParse {
    using namespace boost::units;
    typedef quantity<si::time, double>        time_units_t;
    typedef quantity<si::resistance, double>  resistance_units_t;
    typedef quantity<si::capacitance, double> capacitance_units_t;
    typedef quantity<si::inductance, double>  inductance_units_t;
    struct spef {
      std::string standard;          // i.e. IEEE 1481-1998 - the only one we support for now
      std::string name;
      boost::posix_time::ptime date; // a date *and* time, actually
      time_units_t        t_unit;     // multiplier for all duration values
      capacitance_units_t c_unit;     // capacitance
      resistance_units_t  r_unit;     // and for resistance
      inductance_units_t  l_unit;     // and inductance
    };
  }
}

BOOST_FUSION_ADAPT_STRUCT(
  EDASkel::SpefParse::spef,
  (std::string, standard)
  (std::string, name)
  (EDASkel::SpefParse::time_units_t,        t_unit)
  (EDASkel::SpefParse::capacitance_units_t, c_unit)
  (EDASkel::SpefParse::resistance_units_t,  r_unit)
  (EDASkel::SpefParse::inductance_units_t,  l_unit)
)

#endif // PARSER_SPEFTYPES_H
