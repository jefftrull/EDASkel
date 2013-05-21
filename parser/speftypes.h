#include <boost/fusion/include/adapt_struct.hpp>
#include <string>


#ifndef PARSER_SPEFTYPES_H
#define PARSER_SPEFTYPES_H

namespace EDASkel {
  namespace SpefParse {
    struct spef {
      std::string name;
    };
  }
}

BOOST_FUSION_ADAPT_STRUCT(
  EDASkel::SpefParse::spef,
  (std::string, name)
)

#endif // PARSER_SPEFTYPES_H
