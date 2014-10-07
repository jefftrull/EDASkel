// Basic parser for Standard Parasitic Extraction Format files
// Copyright (C) 2013 Jeffrey Elliot Trull <edaskel@att.net>
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

#ifndef PARSER_SPEF_HPP
#define PARSER_SPEF_HPP

#include "speftypes.h"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

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

    template<typename Iterator, typename SpefVisitor>
    struct spefparser : boost::spirit::qi::grammar<Iterator,
                                                   spef(),
                                                   spefskipper<Iterator> >
    {
      typedef typename SpefVisitor::name_token_value_t name_token_value_t;

      spefparser(SpefVisitor& v) : spefparser::base_type(spef_file), visitor_(v)
      {
        using namespace boost::spirit::qi;
        namespace phx = boost::phoenix;
        using boost::phoenix::val;
        using boost::phoenix::construct;
        using boost::spirit::_1;
        using boost::spirit::_2;
        using boost::spirit::_3;
        using boost::spirit::_4;
        using boost::spirit::_5;
        using boost::spirit::_6;
        using boost::spirit::_a;
        using boost::spirit::_val;

        quoted_string = '"' >> no_skip[*(char_ - '"')] >> '"' ;
        standard    = omit[lexeme["*SPEF"]] >> quoted_string ;
        design_name = omit[lexeme["*DESIGN"]] >> quoted_string ;
        datestr     = omit[lexeme["*DATE"]] >> quoted_string ;
        vendor      = omit[lexeme["*VENDOR"]] >> quoted_string ;
        program     = omit[lexeme["*PROGRAM"]] >> quoted_string ;
        version     = omit[lexeme["*VERSION"]] >> quoted_string ;

        nonspace_str = lexeme[+(char_ - '"' - ' ')];
        design_flow_entry = '"' >> nonspace_str >> (no_skip[" " >> nonspace_str] | attr("")) >> '"' ;
        design_flow = omit[lexeme["*DESIGN_FLOW"]] >> *design_flow_entry ;

        divider     = omit[lexeme["*DIVIDER"]] >> char_ ;
        delimiter   = omit[lexeme["*DELIMITER"]] >> char_ ;
        bus_delimiter = omit[lexeme["*BUS_DELIMITER"]] >> char_ >> char_ ;

        eng_prefixes.add("K", 1E3)("M", 1E-3)("U", 1E-6)("N", 1E-9)("P", 1E-12)("F", 1E-15);
        eng_prefix = lexeme[eng_prefixes] | attr(1.0) ;  // defaults to 1

        t_unit      = omit[lexeme["*T_UNIT"]] >>
                      double_[_a = _1] >>
                      eng_prefix[_val = _a * _1 * val(si::seconds)] >> 'S' ;
        c_unit      = omit[lexeme["*C_UNIT"]] >>
                      double_[_a = _1] >>
                      eng_prefix[_val = _a * _1 * val(si::farads)] >> 'F' ;
        r_unit      = omit[lexeme["*R_UNIT"]] >>
                      double_[_a = _1] >>
                      eng_prefix[_val = _a * _1 * val(si::ohms)] >> "OHM" ;
        l_unit      = omit[lexeme["*L_UNIT"]] >>
                      double_[_a = _1] >>
                      eng_prefix[_val = _a * _1 * val(si::henrys)] >> (lit("HENRY") | 'H') ;

        netname = lexeme[+(char_("a-zA-Z0-9[]/:_") | ('\\' >> char_("[]/")))];
        name_map_entry = ('*' >> lexeme[+ascii::digit] >> netname)[
           phx::bind(name_map_symtab.add, _1,
                     phx::bind(&SpefVisitor::name_map_entry, phx::ref(visitor_), _2))] ;
        name_map = omit[lexeme["*NAME_MAP"]] >> *name_map_entry ;

        port_def = (lexeme['*' >> name_map_symtab] >> char_("IOB") >>
                    // stuff I don't understand yet
                    *omit[("*C" >> lexeme[double_] >> lexeme[double_]) |
                          ("*L" >> double_) |
                          ("*S" >> double_ >> double_)])[
                            phx::bind(&SpefVisitor::port_definition, phx::ref(visitor_), _1, _2)] ;
        ports = lexeme["*PORTS"] >> *port_def;

        connection = '*' >> ((lit('P') >> '*' >> name_map_symtab)[
                               phx::bind(&SpefVisitor::net_port_connection, phx::ref(visitor_), _r1, _1)] |
                             (lit('I') >> '*' >> name_map_symtab >> ':' >> as_string[lexeme[+ascii::alnum]])[
                               phx::bind(&SpefVisitor::net_inst_connection, phx::ref(visitor_), _r1, _1, _2)])
                         >> char_("IOB")
                         >> -(lit("*C") >> double_ >> double_)   // location in X/Y coordinates
                         >> -(lit("*L") >> double_)              // pin load
                         >> -(lit("*D") >> lexeme[+(ascii::alnum | '_')]) ;

        gndcapline = (uint_
                      >> '*' >> name_map_symtab >> as_string[(':' >> lexeme[+ascii::alnum]) | attr("")]
                      >> double_)[
                        phx::bind(&SpefVisitor::cgnd, phx::ref(visitor_),
                                  _r1, _1, _2, _3, _4 * phx::cref(c_unit_value))];

        capline = (uint_
                   >> '*' >> name_map_symtab >> as_string[(':' >> lexeme[+ascii::alnum]) | attr("")]
                   >> '*' >> name_map_symtab >> as_string[(':' >> lexeme[+ascii::alnum]) | attr("")]
                   >> double_)[
                     phx::bind(&SpefVisitor::capacitor, phx::ref(visitor_),
                               _r1, _1, _2, _3, _4, _5, _6 * phx::cref(c_unit_value))];

        resline = (uint_ >> '*' >> name_map_symtab >> as_string[(':' >> lexeme[+ascii::alnum]) | attr("")]
                         >> '*' >> name_map_symtab >> as_string[(':' >> lexeme[+ascii::alnum]) | attr("")]
                         >> double_)[
                           phx::bind(&SpefVisitor::resistor, phx::ref(visitor_),
                                     _r1, _1, _2, _3, _4, _5, _6 * phx::cref(r_unit_value))];

        net_def =  lit("*D_NET") >>
                   ('*' >> name_map_symtab >> double_)[
                     _a = _1,
                     phx::bind(&SpefVisitor::net_definition, phx::ref(visitor_),
                               _1, _2 * phx::cref(c_unit_value))] >>
                   -("*CONN" >> *connection(_a)) >>
                   -("*CAP" >> *(capline(_a) | gndcapline(_a))) >>
                   -("*RES" >> *resline(_a)) >>
                   "*END" ;

        nets = *net_def ;

        spef_file %= -omit[lexeme["SPEF"]] >> standard >> design_name   // TODO: any order?
                                           >> omit[datestr] >> vendor >> program >> version
                                           >> design_flow
                                           >> omit[divider] >> omit[delimiter] >> omit[bus_delimiter]
                                           >> t_unit
                                           >> omit[c_unit[phx::ref(c_unit_value) = _1]]
                                           >> omit[r_unit[phx::ref(r_unit_value) = _1]]
                                           >> l_unit
                                           >> -name_map >> -ports >> nets;

        spef_file.name("SPEF top level");
        standard.name("SPEF standard version");
        design_name.name("DESIGN name");
        datestr.name("Date");
        vendor.name("Vendor");
        program.name("Program");
        version.name("Version");
        design_flow.name("Design Flow");

        t_unit.name("time unit declaration");
        c_unit.name("capacitance unit declaration");
        r_unit.name("resistance unit declaration");
        l_unit.name("inductance unit declaration");

        port_def.name("Port Definition");
        ports.name("Port List");
        net_def.name("Net Definition");
        nets.name("Net List");

        connection.name("net connection definition");
        capline.name("capacitor definition");
        gndcapline.name("grounded capacitor definition");
        resline.name("resistor definition");

        BOOST_SPIRIT_DEBUG_NODE(spef_file);
        BOOST_SPIRIT_DEBUG_NODE(standard);
        BOOST_SPIRIT_DEBUG_NODE(design_name);
        BOOST_SPIRIT_DEBUG_NODE(datestr);
        BOOST_SPIRIT_DEBUG_NODE(vendor);
        BOOST_SPIRIT_DEBUG_NODE(program);
        BOOST_SPIRIT_DEBUG_NODE(version);
        BOOST_SPIRIT_DEBUG_NODE(design_flow);
        BOOST_SPIRIT_DEBUG_NODE(design_flow_entry);
        BOOST_SPIRIT_DEBUG_NODE(nonspace_str);
        BOOST_SPIRIT_DEBUG_NODE(eng_prefix);
        BOOST_SPIRIT_DEBUG_NODE(t_unit);
        BOOST_SPIRIT_DEBUG_NODE(r_unit);
        BOOST_SPIRIT_DEBUG_NODE(c_unit);
        BOOST_SPIRIT_DEBUG_NODE(name_map_entry);
        BOOST_SPIRIT_DEBUG_NODE(name_map);
        BOOST_SPIRIT_DEBUG_NODE(port_def);
        BOOST_SPIRIT_DEBUG_NODE(ports);
        BOOST_SPIRIT_DEBUG_NODE(net_def);
        BOOST_SPIRIT_DEBUG_NODE(nets);
        BOOST_SPIRIT_DEBUG_NODE(connection);
        BOOST_SPIRIT_DEBUG_NODE(capline);
        BOOST_SPIRIT_DEBUG_NODE(gndcapline);
        BOOST_SPIRIT_DEBUG_NODE(resline);

        on_error<fail>(spef_file, std::cerr << val("Error! Expecting ")
                                            << boost::spirit::_4
                                            << val(" here: \"")
                                            << construct<std::string>(boost::spirit::_3,
                                                                      boost::spirit::_2)
                                            << val("\"") << std::endl);

      }
      
    private:

      typedef spefskipper<Iterator> skipper_t;
      template<typename Signature>
      using Rule = boost::spirit::qi::rule<Iterator, Signature, skipper_t>;
      typedef boost::spirit::qi::rule<Iterator, skipper_t> no_attr_rule_t;

      Rule<spef()> spef_file;

      Rule<std::string()> quoted_string;
      Rule<std::string()> design_name, standard, datestr, vendor, program, version;

      typename boost::spirit::qi::rule<Iterator, std::string()> nonspace_str;
      typedef std::pair<std::string, std::string> design_flow_entry_t;
      Rule<design_flow_entry_t()> design_flow_entry;
      Rule<design_flow_map_t()> design_flow;

      boost::spirit::qi::symbols<char, double> eng_prefixes;
      Rule<double()> eng_prefix;

      Rule<char()> divider, delimiter;
      Rule<std::pair<char, char> > bus_delimiter;

      template <typename Quantity>
      using unit_rule = typename boost::spirit::qi::rule<Iterator,
                                                         quantity<Quantity, double>(),
                                                         skipper_t,
                                                         boost::spirit::locals<double> >;
      unit_rule<si::time>                t_unit;
      unit_rule<si::resistance>          r_unit;
      quantity<si::resistance, double>   r_unit_value;
      unit_rule<si::capacitance>         c_unit;
      quantity<si::capacitance, double>  c_unit_value;
      unit_rule<si::inductance>          l_unit;

      typename boost::spirit::qi::rule<Iterator, std::string()> netname;
      no_attr_rule_t name_map_entry, name_map;   // semantic actions create symtab entries
      boost::spirit::qi::symbols<char, name_token_value_t> name_map_symtab;

      no_attr_rule_t port_def, ports;

      typedef boost::spirit::qi::rule<Iterator, skipper_t, void(name_token_value_t)>
          no_attr_net_rule_t;
      no_attr_net_rule_t connection, capline, gndcapline, resline;
      boost::spirit::qi::rule<Iterator, skipper_t,
                              boost::spirit::locals<name_token_value_t> > net_def;
      no_attr_rule_t nets;

      SpefVisitor & visitor_;      // user-defined object that responds to data

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
