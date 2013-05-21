// Basic parser for Standard Parasitic Extraction Format files

#ifndef PARSER_SPEF_HPP
#define PARSER_SPEF_HPP

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include "speftypes.h"

namespace EDASkel {
  namespace SpefParse {

    // typedef for stream iterator we will use
    typedef boost::spirit::istream_iterator SpefIter;

    // create skipper for comments and whitespace
    // adapted from a presentation at Boostcon 2010 by Michael Caisse
    template <typename Iterator>
    struct spefskipper : boost::spirit::qi::grammar< Iterator >
    {
      spefskipper() : spefskipper::base_type(skip_it)
      {
        using namespace boost::spirit::qi;

        comment = "//" >> *( char_ - eol ) >> eol ;
        skip_it = comment | space ;

      }
      boost::spirit::qi::rule<Iterator> skip_it;
      boost::spirit::qi::rule<Iterator> comment;
    };

    // storage for node number aliases
    typedef std::map<std::size_t, std::string> alias_map_t;

    template<typename Iterator>
    struct spefparser : boost::spirit::qi::grammar<Iterator,
                                                   spef(),
                                                   spefskipper<Iterator> >
    {
      spefparser() : spefparser::base_type(spef_file)
      {
        using namespace boost::spirit::qi;
        using boost::phoenix::val;
        using boost::phoenix::construct;

        design_name = omit[lexeme["*DESIGN"]] >> '"' >> no_skip[*(char_ - '"')] >> '"' ;

        spef_file = design_name ;                  // the only thing we parse for now

        spef_file.name("SPEF top level");
        design_name.name("DESIGN name");

        on_error<fail>(spef_file, std::cerr << val("Error! Expecting ")
                                            << boost::spirit::_4
                                            << val(" here: \"")
                                            << construct<std::string>(boost::spirit::_3,
                                                                      boost::spirit::_2)
                                            << val("\"") << std::endl);

      }
      
      typedef spefskipper<Iterator> skipper_t;
      template<typename Signature>
      struct Rule
      {
        typedef boost::spirit::qi::rule<Iterator, Signature, skipper_t> type;
      };

      typename Rule<spef()>::type spef_file;
      typename Rule<std::string()>::type design_name;
    };
  }  // namespace SpefParse
} // namespace EDASkel


// grammar taken from Wikipedia

// from the header, take design name, extraction tool, and units
// *R_UNIT, *C_UNIT but there is also *T_UNIT and *L_UNIT
// also *DIVIDER and *DELIMITER values will affect the parsing
// non alnums are supposed to be escaped with the exception of the divider and delimiter

// name map - stick results into symbol table

// *NAME_MAP
// *509 F_C_EP2
// *513 TOP/BUF_ZCLK_2_pin_Z_1

// now can use F_C_EP2 or *509 wherever a node name is required

// *PORTS
// *1 I
// *2 I
// *3 O
// *4 B   << I,O, or B for in, out, bidir

// *D_NET regcontrol_top/GRC/n13345 1.94482   << lumped C
// *CONN
// *I regcontrol_top/GRC/U9743:E I   << can also be *P for top-level port
// *C 537.855 9150.11     << XY coordinate
// *L 3.7000   << pin cap - if output, this will be *D and driving cell
// *CAP
// 1 regcontrol_top/GRC/U9743:E 0.936057   << lumped to ground style
// 2 regcontrol_top/GRC/U9709:A regcontrol_top/GRC/U10716:Z 0.622675  << coupling cap style
// *RES
// 1 regcontrol_top/GRC/U9743:E regcontrol_top/GRC/U9407:Z 10.7916
// *END  << of this particular D_NET

// semantic checks for SPEF
// lumped cap vs. total of parasitics
// loop detection (and possibly removal) via graph
// compare netlists vs. a design - a good place for applying generic programming
// (spef as netlist, DEF as netlist, vlog as netlist...)

// standard expects 1 (nominal or subject to interpretation) or 3 (min:typ:max) values for caps/resistors
// but in the future may have any number of corners separated by colons

#endif // PARSER_SPEF_HPP
