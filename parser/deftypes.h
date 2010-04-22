// class definitions for DEF parsing results
// Copyright 2010 Jeffrey Elliot Trull <linmodemstudent@gmail.com>
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
