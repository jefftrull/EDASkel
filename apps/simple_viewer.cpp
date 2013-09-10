// A simple DEF design viewer, part of EDASkel, a sample EDA app
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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <QApplication>
#include <QGraphicsView>

#include <boost/program_options.hpp>

#include "../parser/lefparser.h"
#include "../parser/lefsem.h"
#include "../parser/defparser.h"
#include "../parser/defsem.h"
#include "../db/simpledb.h"
#include "../gui/designscene.h"
#include "../gui/designview.h"

namespace EDASkel {
  extern DefParse::DefTokens<LefDefLexer> defTokens;
  extern DefParse::defparser<DefParse::DefTokens<LefDefLexer>::iterator_type,
			     DefParse::DefTokens<LefDefLexer>::lexer_def > defParser;
  extern LefParse::LefTokens<LefDefLexer> lefTokens;
  extern LefParse::lefparser<LefParse::LefTokens<LefDefLexer>::iterator_type,
			     LefParse::LefTokens<LefDefLexer>::lexer_def > lefParser;
}

int main(int argc, char **argv) {

  using namespace LefParse;
  using namespace DefParse;
  using namespace SimpleDB;
  using namespace EDASkel;
  using namespace boost::spirit::qi;
  using boost::spirit::qi::space;

  // parse command line options

  QApplication app(argc, argv);       // let Qt take any args it recognizes

  std::string lef_fn, def_fn;         // input filenames

  namespace po = boost::program_options;
  try {
     po::options_description desc("Simple DEF Viewer options");
     desc.add_options()
        ( "help,h",                                              "this help message" )
        ( "lef",    po::value<std::string>(&lef_fn)->required(), "LEF file to use" )
        ( "def",    po::value<std::string>(&def_fn)->required(), "DEF to read and display" ) ;

     po::variables_map vm;
     po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

     if (vm.count("help")) {
        std::cout << "usage: " << argv[0] << " --lef filename.lef --def filename.def" << std::endl;
        return 0;
     }

     po::notify(vm);
  } catch (std::exception const& e) {
     std::cerr << "Argument Error: " << e.what() << std::endl;
     return 1;
  }

  // read LEF/DEF files into library and database
  std::ifstream lefin(lef_fn.c_str());
  if (!lefin.is_open()) {
    std::cerr << "Failed to open input LEF file " << lef_fn << std::endl;
    return 1;
  }
  // do not skip whitespace
  lefin.unsetf(std::ios::skipws);
  LefDefIter beg = LefDefIter(lefin), end;
  LefTokens<LefDefLexer>::iterator_type lef_it = lefTokens.begin(beg, LefDefIter());
  LefTokens<LefDefLexer>::iterator_type lef_end = lefTokens.end();
  lef lef_ast;
  if (!parse(lef_it, lef_end, lefParser, lef_ast) ||
      (lef_it != lef_end)) {
    std::cerr << "LEF parse failed\n";
    if (lef_it != lef_end) {
       std::cerr << "did not consume all input; ";
       std::cerr << std::distance(lef_it, lef_end) << " tokens remain:" << std::endl;
       std::transform(lef_it, lef_end,
                      std::ostream_iterator<boost::iterator_range<LefDefIter> >(std::cerr, "|"),
                      [](LefTokens<LefDefLexer>::iterator_type::value_type const& tok) {
                        // token's value is a variant which initially stores the original iterator
                        // range matched by the token, i.e., a range of characters
                           return boost::get<boost::iterator_range<LefDefIter> >(tok.value());
                      });
    }
    return 1;
  }
  LefChecker<Library> lefchk;
  Library lib;
  if (!lefchk.CheckAndInsert(lef_ast, lib)) {
    std::cerr << "LEF semantic check failed\n";
    return 1;
  }

  std::ifstream defin(def_fn.c_str());
  if (!defin.is_open()) {
    std::cerr << "Failed to open input DEF file " << def_fn << std::endl;
    return 1;
  }
  defin.unsetf(std::ios::skipws);
  beg = LefDefIter(defin);
  DefTokens<LefDefLexer>::iterator_type it = defTokens.begin(beg, end);
  DefTokens<LefDefLexer>::iterator_type lex_end = defTokens.end();
  def def_ast;
  if (!parse(it, lex_end, defParser, def_ast) ||
      (beg != end)) {
    std::cerr << "DEF parse failed\n";
    if (it != lex_end) {
       std::cerr << "did not consume all input; ";
       std::cerr << std::distance(it, lex_end) << " tokens remain:" << std::endl;
       std::transform(it, lex_end,
                      std::ostream_iterator<boost::iterator_range<LefDefIter> >(std::cerr, "|"),
                      [](DefTokens<LefDefLexer>::iterator_type::value_type const& tok) {
                           return boost::get<boost::iterator_range<LefDefIter> >(tok.value());
                      });
    }
    return 1;
  }
  DefChecker<Database, Library> defchk;
  Database db;
  if (!defchk.CheckAndInsert(def_ast, lib, db)) {
    std::cerr << "DEF semantic checks failed\n";
    return 1;
  }

  boost::iterator_range<Database::InstIter> instRange = db.getInstances();
  std::cerr << "DEF contains " << std::distance(std::begin(instRange), std::end(instRange)) << " instances\n";

  // create a "scene" (model) for display from the filled-out library and database
  DesignScene<Database, Library> myScene(db, lib);

  // Our special view adds some keyboard shortcuts and shows the design
  // "right side up"
  DesignView view(&myScene);
  view.fitInView(myScene.sceneRect(), Qt::KeepAspectRatio);
  view.show();
  app.exec();

  return 0;

}
