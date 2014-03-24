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

#include "../tclint/qttclnotifier.h"

namespace EDASkel {
  extern DefParse::DefTokens<LefDefLexer> defTokens;
  extern DefParse::defparser<DefParse::DefTokens<LefDefLexer>::iterator_type,
			     DefParse::DefTokens<LefDefLexer>::lexer_def > defParser;
  extern LefParse::LefTokens<LefDefLexer> lefTokens;
  extern LefParse::lefparser<LefParse::LefTokens<LefDefLexer>::iterator_type,
			     LefParse::LefTokens<LefDefLexer>::lexer_def > lefParser;
}

struct SetupException : std::exception {
  explicit SetupException(std::string err) : err_(err) {}
  const char* what() const noexcept { return err_.c_str(); }
  virtual ~SetupException() noexcept {}
private:
  std::string err_;
};

// Container for the Tcl commands we define, and any data they require
// Handles cleanup also
struct CommandState {
  CommandState(Tcl_Interp* interp) : interp_(interp) {
    // get arc and argv from interp
    // (Tcl_Main does not supply them to us, but does create Tcl variables of the same name)
    Tcl_Obj* argv_obj = Tcl_GetVar2Ex(interp_, "argv", NULL, TCL_GLOBAL_ONLY);
    if (argv_obj == NULL) {
      throw SetupException("could not get value of argv");
    }
    Tcl_Obj** objv;
    int argc;
    if (Tcl_ListObjGetElements(interp_, argv_obj, &argc, &objv) != TCL_OK) {
      throw SetupException("could not extract argv");
    }

    // now objv is an array of Tcl_Obj*.  We must turn each of those into strings
    using namespace std;
    vector<string> argv;
    transform(objv, objv+argc, back_inserter(argv), Tcl_GetString);
      
    // parse command line options
    string lef_fn, def_fn;         // input filenames

    // no try block here - let failures propagate out of ctor
    namespace po = boost::program_options;
    po::options_description desc("Simple DEF Viewer options");
    desc.add_options()
      ( "help,h",                                              "this help message" )
      ( "lef",    po::value<string>(&lef_fn)->required(), "LEF file to use" )
      ( "def",    po::value<string>(&def_fn)->required(), "DEF to read and display" ) ;

    po::variables_map vm;
    po::store(po::command_line_parser(argv).options(desc).run(), vm);

    if (vm.count("help")) {
      throw SetupException("usage: " + argv[0] + " --lef filename.lef --def filename.def");
    }

    po::notify(vm);

    // TODO consume our args from Tcl's variables?

    using namespace LefParse;
    using namespace DefParse;
    using namespace SimpleDB;
    using namespace EDASkel;
    using namespace boost::spirit::qi;
    using boost::spirit::qi::space;

    // read LEF/DEF files into library and database
    ifstream lefin(lef_fn.c_str());
    if (!lefin.is_open()) {
      throw SetupException("Failed to open input LEF file " + lef_fn);
    }
    // do not skip whitespace
    lefin.unsetf(ios::skipws);
    LefDefIter beg = LefDefIter(lefin), end;
    LefTokens<LefDefLexer>::iterator_type lef_it = lefTokens.begin(beg, LefDefIter());
    LefTokens<LefDefLexer>::iterator_type lef_end = lefTokens.end();
    lef lef_ast;
    if (!parse(lef_it, lef_end, lefParser, lef_ast) ||
        (lef_it != lef_end)) {
      if (lef_it != lef_end) {
        cerr << "did not consume all input; ";
        cerr << distance(lef_it, lef_end) << " tokens remain:" << endl;
        transform(lef_it, lef_end,
                       ostream_iterator<boost::iterator_range<LefDefIter> >(cerr, "|"),
                       [](LefTokens<LefDefLexer>::iterator_type::value_type const& tok) {
                         // token's value is a variant which initially stores the original iterator
                         // range matched by the token, i.e., a range of characters
                         return boost::get<boost::iterator_range<LefDefIter> >(tok.value());
                       });
      }
      throw SetupException("LEF parse failed");
    }
    LefChecker<Library> lefchk;
    if (!lefchk.CheckAndInsert(lef_ast, lib_)) {
      throw SetupException("LEF semantic check failed");
    }

    ifstream defin(def_fn.c_str());
    if (!defin.is_open()) {
      throw SetupException("Failed to open input DEF file " + def_fn);
    }
    defin.unsetf(ios::skipws);
    beg = LefDefIter(defin);
    DefTokens<LefDefLexer>::iterator_type it = defTokens.begin(beg, end);
    DefTokens<LefDefLexer>::iterator_type lex_end = defTokens.end();
    def def_ast;
    if (!parse(it, lex_end, defParser, def_ast) ||
        (beg != end)) {
      if (it != lex_end) {
        cerr << "did not consume all input; ";
        cerr << distance(it, lex_end) << " tokens remain:" << endl;
        transform(it, lex_end,
                       ostream_iterator<boost::iterator_range<LefDefIter> >(cerr, "|"),
                       [](DefTokens<LefDefLexer>::iterator_type::value_type const& tok) {
                         return boost::get<boost::iterator_range<LefDefIter> >(tok.value());
                       });
      }
      throw SetupException("DEF parse failed");
    }
    DefChecker<Database, Library> defchk;
    if (!defchk.CheckAndInsert(def_ast, lib_, db_)) {
      throw SetupException("DEF semantic checks failed");
    }

    boost::iterator_range<Database::InstIter> instRange = db_.getInstances();
    cerr << "DEF contains " << distance(std::begin(instRange), std::end(instRange)) << " instances\n";

    // create a "scene" (model) for display from the filled-out library and database
    scene_ = new DesignScene(db_, lib_);
    auto design_extent = scene_->sceneRect();

    // Our special view adds some keyboard shortcuts and shows the design
    // "right side up"
    DesignView* view = new DesignView(scene_);  // owner??
    view->fitInView(design_extent, Qt::KeepAspectRatio);
    view->show();

    // register static methods using "clientData" to hold the "this" pointer
    Tcl_CreateObjCommand(interp_, "instances", &CommandState::instances, this, NULL);
    Tcl_CreateObjCommand(interp_, "highlight", &CommandState::highlight, this, NULL);

    // Ensure final cleanup invokes our destructor
    Tcl_CreateExitHandler(&CommandState::cleanup, this);
  }

private:
  // these statics provide a C-style function pointer interface that Tcl expects,
  // but redirect to regular instance methods via the ClientData pointer
  static int instances(ClientData me, Tcl_Interp*, int objc, Tcl_Obj* const objv[]) {
    // strategy TODO: should user syntax/semantic checking be in Tcl, here, or in the instance method?
    return static_cast<CommandState*>(me)->instances(objc, objv);
  }
  static int highlight(ClientData me, Tcl_Interp*, int objc, Tcl_Obj* const objv[]) {
    return static_cast<CommandState*>(me)->highlight(objc, objv);
  }

  // arranges for destructor to be called when interpreter is shut down
  static void cleanup(ClientData me) {
    delete static_cast<CommandState*>(me);
  }

  // implementations
  int instances(int objc, Tcl_Obj* const objv[]) {
    // syntax check here for now
    if (objc != 1) {
      Tcl_WrongNumArgs(interp_, 1, objv, NULL);
      return TCL_ERROR;
    }

    std::vector<Tcl_Obj*> result_objv;
    auto instRange = db_.getInstances();
    std::transform(std::begin(instRange), std::end(instRange), back_inserter(result_objv),
                   [](std::shared_ptr<SimpleDB::Instance> i) {
                     return Tcl_NewStringObj(i->getName().c_str(), i->getName().length());
                   });
    Tcl_SetObjResult(interp_, Tcl_NewListObj(result_objv.size(), result_objv.data()));

    return TCL_OK;
  }

  int highlight(int objc, Tcl_Obj* const objv[]) {
    if (objc != 2) {
      Tcl_WrongNumArgs(interp_, 1, objv, "instance_name");
      return TCL_ERROR;
    }

    // attempt to highlight
    if (!scene_->highlightInstance(Tcl_GetString(objv[1]))) {
      Tcl_AppendResult(interp_, "could not find instance ", Tcl_GetString(objv[1]), NULL);
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  Tcl_Interp* interp_;
  SimpleDB::Database db_;
  SimpleDB::Library lib_;
  DesignScene* scene_;

};

// provide former contents of main() in a function with Tcl_AppInit signature 
int startup(Tcl_Interp* interp) {
  // get argc and argv back from interp

  // instantiate "notifier" to combine Tcl and Qt events
  QtTclNotify::QtTclNotifier::setup();

  // define commands
  try {
    new CommandState(interp);   // NOT a leak. Tcl interp owns data
  } catch(std::exception& e) {
    std::cerr << e.what() << std::endl;
    Tcl_Exit(1);
  }

  // run Qt's event loop
  Tcl_SetMainLoop([]() { QApplication::exec(); });

  return TCL_OK;
}


int main(int argc, char **argv) {

  QApplication app(argc, argv);       // let Qt take any args it recognizes

  // create a Tcl interpreter, connect it to the terminal, and launch our code
  Tcl_Main(argc, argv, startup);

}
