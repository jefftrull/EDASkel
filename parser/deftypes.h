// class definitions for DEF parsing results
// Copyright 2010 Jeffrey Elliot Trull <linmodemstudent@gmail.com>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>
#include <string>

struct defplcinfo {
  std::string plcfix;       // BOZO make boolean "fixed" or something
  int x, y;
  std::string orient;       // BOZO make this enum
};

BOOST_FUSION_ADAPT_STRUCT(
  defplcinfo,
  (std::string, plcfix)
  (int, x)
  (int, y)
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

