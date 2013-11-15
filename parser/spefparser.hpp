// Basic parser for Standard Parasitic Extraction Format files

#ifndef PARSER_SPEF_HPP
#define PARSER_SPEF_HPP

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include "speftypes.h"

namespace EDASkel {
  namespace SpefParse {

    using namespace boost::units;

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
        using boost::spirit::_1;
        using boost::spirit::_a;
        using boost::spirit::_val;

        eng_prefixes.add("K", 1E3)("M", 1E-3)("U", 1E-6)("N", 1E-9)("P", 1E-12)("F", 1E-15);
        eng_prefix = lexeme[eng_prefixes] | (eps >> attr(1.0)) ;  // defaults to 1

        design_name = omit[lexeme["*DESIGN"]] >> '"' >> no_skip[*(char_ - '"')] >> '"' ;
        standard    = omit[lexeme["*SPEF"]] >> '"' >> no_skip[*(char_ - '"')] >> '"' ;
        t_unit      = omit[lexeme["*T_UNIT"]] >>
                      double_[_a = _1] >>
                      eng_prefix[_val = _a * _1 * val(si::seconds)] >> 'S' ;
        r_unit      = omit[lexeme["*R_UNIT"]] >>
                      double_[_a = _1] >>
                      eng_prefix[_val = _a * _1 * val(si::ohms)] >> "OHM" ;
        c_unit      = omit[lexeme["*C_UNIT"]] >>
                      double_[_a = _1] >>
                      eng_prefix[_val = _a * _1 * val(si::farads)] >> 'F' ;

        spef_file = omit[lexeme["SPEF"]] >> standard >> design_name   // TODO: any order?
                                         >> t_unit >> r_unit >> c_unit;

        spef_file.name("SPEF top level");
        design_name.name("DESIGN name");
        standard.name("SPEF standard version");
        t_unit.name("time unit declaration");
        r_unit.name("resistance unit declaration");
        c_unit.name("capacitance unit declaration");

        BOOST_SPIRIT_DEBUG_NODE(spef_file);
        BOOST_SPIRIT_DEBUG_NODE(design_name);
        BOOST_SPIRIT_DEBUG_NODE(standard);
        BOOST_SPIRIT_DEBUG_NODE(eng_prefix);
        BOOST_SPIRIT_DEBUG_NODE(t_unit);
        BOOST_SPIRIT_DEBUG_NODE(r_unit);
        BOOST_SPIRIT_DEBUG_NODE(c_unit);

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
      typename Rule<std::string()>::type design_name, standard;
      boost::spirit::qi::symbols<char, double> eng_prefixes;
      typename Rule<double()>::type eng_prefix;

      template <typename Quantity>
      struct unit_rule {
        typedef typename boost::spirit::qi::rule<Iterator, quantity<Quantity, double>(), skipper_t,
                                        boost::spirit::locals<double> > type;
      };
      typename unit_rule<si::time>::type        t_unit;
      typename unit_rule<si::resistance>::type  r_unit;
      typename unit_rule<si::capacitance>::type c_unit;
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
