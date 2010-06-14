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

#include <boost/spirit/include/support_multi_pass.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tr1/tuple>

#include <QApplication>
#include <QGraphicsView>

#include "../parser/lefparser.h"
#include "../parser/lefsem.h"
#include "../parser/defparser.h"
#include "../parser/defsem.h"
#include "../db/simpledb.h"
#include "../gui/designscene.h"
#include "../gui/designview.h"

typedef std::istreambuf_iterator<char> iterator_base;
typedef boost::spirit::multi_pass<iterator_base> mp_stream_iterator;

namespace EDASkel {
  extern DefParse::defparser<mp_stream_iterator> defFileParser;
  extern LefParse::lefparser<mp_stream_iterator> lefFileParser;
  extern lefdefskipper<mp_stream_iterator> lefdefFileSkipper;
}

int main(int argc, char **argv) {

  using namespace LefParse;
  using namespace DefParse;
  using namespace SimpleDB;
  using namespace EDASkel;
  using namespace boost::spirit::qi;
  using boost::spirit::qi::space;

  // parse command line options
  // I want to do this cleanly, but can't decide between getopt (kind of ugly,
  // not cross-platform) and Boost program_options (non header-only).
  // doing it the ugly manual way...

  QApplication app(argc, argv);       // let Qt take any args it recognizes
  QStringList args = app.arguments();
  boost::optional<std::string> lef_fn, def_fn;
  for (int i = 0; i < (args.size() - 1); ++i) {
    std::string arg(args.at(i).toLocal8Bit().constData());
    if (arg == "--lef")
      lef_fn = args.at(++i).toLocal8Bit().constData();
    else if (arg == "--def")
      def_fn = args.at(++i).toLocal8Bit().constData();
  }
  if (!lef_fn) {
    std::cerr << "Please specify a LEF file with the --lef option\n";
    return 1;
  }
  if (!def_fn) {
    std::cerr << "Please specify a DEF file with the --def option\n";
    return 1;
  }

  // read LEF/DEF files into library and database
  std::ifstream lefin(lef_fn->c_str());
  if (!lefin.is_open()) {
    std::cerr << "Failed to open input LEF file " << *lef_fn << std::endl;
    return 1;
  }
  mp_stream_iterator beg = boost::spirit::make_default_multi_pass(iterator_base(lefin)), end;
  lef lef_ast;
  if (!phrase_parse(beg, end, lefFileParser, lefdefFileSkipper, lef_ast) ||
      (beg != end)) {
    std::cerr << "LEF parse failed\n";
    if (beg != end)
      std::cerr << "did not consume all input; remaining is:\n" << std::string(beg, end) << std::endl;
    return 1;
  }
  LefChecker<Library> lefchk;
  Library lib;
  if (!lefchk.CheckAndInsert(lef_ast, lib)) {
    std::cerr << "LEF semantic check failed\n";
    return 1;
  }

  std::ifstream defin(def_fn->c_str());
  if (!defin.is_open()) {
    std::cerr << "Failed to open input DEF file " << *def_fn << std::endl;
    return 1;
  }
  beg = mp_stream_iterator(defin);
  def def_ast;
  if (!phrase_parse(beg, end, defFileParser, lefdefFileSkipper, def_ast) ||
      (beg != end)) {
    std::cerr << "DEF parse failed\n";
    if (beg != end)
      std::cerr << "remaining input is as follows:\n" << std::string(beg, end) << std::endl;
    return 1;
  }
  DefChecker<Database, Library> defchk;
  Database db;
  if (!defchk.CheckAndInsert(def_ast, lib, db)) {
    std::cerr << "DEF semantic checks failed\n";
    return 1;
  }

  typedef Database::InstIter InstIter;
  InstIter iit, instend;
  std::tr1::tie(iit, instend) = db.getInstances();
  std::cerr << "DEF contains " << distance(iit, instend) << " instances\n";

  // create a "scene" (model) for display from the filled-out library and database
  DesignScene<Database, Library> myScene(db, lib);
  QRectF sr = myScene.sceneRect();

  // Our special view adds some keyboard shortcuts and shows the design
  // "right side up"
  DesignView view(&myScene);
  view.fitInView(myScene.sceneRect(), Qt::KeepAspectRatio);
  view.show();
  app.exec();

  return 0;

}
