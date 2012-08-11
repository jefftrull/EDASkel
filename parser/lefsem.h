// Semantic analysis for a Boost Spirit-based LEF parser, part of EDASkel, a sample EDA app
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

#include "leftypes.h"

#if !defined(EDASKEL_LEF_SEMANTIC)
#define EDASKEL_LEF_SEMANTIC

namespace LefParse {

  // For a discussion of the error handling strategy I'm using here see "defsem.h"

  // list of known errors
  enum LefCheckError { 
    LEFERR_SITE_REDEF,
    LEFERR_MACRO_REDEF,
    LEFERR_MACRO_NOCLASS
  };

  // default handling approach
  template<LefCheckError>
  class LefCheckPolicy {
  public:
    static bool const silent = false;  // do not be silent (give an error message)
    static bool const skip = true;     // do skip (i.e., do not take action on this element)
    static bool const abort = false;   // do not abort (i.e., try to keep going)
  };

  // and now a member function to check the results of a LEF parse and conditionally insert
  // the elements into a database
  // Look at simpledb.h for examples of the "Concept" Lib must implement
  template<class Lib, template<LefCheckError> class CheckPolicy = LefCheckPolicy>
    class LefChecker {
  public:
  bool CheckAndInsert(const lef& lefin, Lib& lib) {
    // get site definitions from AST, check, and install into lib
    typedef typename Lib::Site LibSite;
    typedef typename Lib::SitePtr LibSitePtr;
    for (std::vector<Site>::const_iterator sit = lefin.sites.begin();
	 sit != lefin.sites.end(); ++sit) {
      // All I can think of to check for sites is whether they've been defined already
      LibSitePtr prevsite = lib.findSite(sit->name);
      bool site_was_redefined = (prevsite != LibSitePtr());
      if (site_was_redefined && !CheckPolicy<LEFERR_SITE_REDEF>::silent) {
	// really need some more diagnostics here: file and line number of original and new definitions
	std::cerr << "LEF Checker: Site " << sit->name << " is being redefined\n";
      }
      if (!site_was_redefined || !CheckPolicy<LEFERR_SITE_REDEF>::skip) {
	// "don't skip" in this case will mean "overwrite previous definition"
	// I don't think having both in the library, with the same name, is a good idea.
	// remove the original one first
	if (site_was_redefined)
	  lib.removeSite(prevsite);

	LibSitePtr sptr(new LibSite(sit->name));
	sptr->setClass(sit->class_);
	if (sit->symmetry)
	  sptr->setSymmetry(*(sit->symmetry));
	sptr->setDimensions(sit->width, sit->height);
	lib.addSite(sptr);
      }
      if (site_was_redefined && CheckPolicy<LEFERR_SITE_REDEF>::abort)
	return false;
    }
    // get cell definitions from AST
    typedef typename Lib::Cell LibCell;
    typedef typename Lib::CellPtr LibCellPtr;
    for (std::vector<lefmacro>::const_iterator mit = lefin.macros.begin();
	 mit != lefin.macros.end(); ++mit) {
      // check to see if this cell is being redefined
      LibCellPtr prevcell = lib.findCell(mit->name);
      bool cell_was_redefined = (prevcell != LibCellPtr());
      if (cell_was_redefined && !CheckPolicy<LEFERR_MACRO_REDEF>::silent) {
	std::cerr << "LEF Checker: Cell " << mit->name << " is being redefined\n";
      }
      if (!cell_was_redefined || !CheckPolicy<LEFERR_MACRO_REDEF>::skip) {
	if (cell_was_redefined)
	  lib.removeCell(prevcell);

	LibCellPtr cptr(new LibCell(mit->name));
	if (mit->class_)
	  cptr->setClass(*(mit->class_));
	else {
	  cptr->setClass(SITECLASS_CORE);
	  if (!CheckPolicy<LEFERR_MACRO_NOCLASS>::silent)
	    // LEF/DEF 5.6 specifies a warning here
	    std::cerr << "LEF Checker: Cell " << mit->name << " has no class defined; assuming CORE\n";
	}
	if (mit->origin)
	  cptr->setOrigin(std::make_pair(mit->origin->x, mit->origin->y));
	else
	  cptr->setOrigin(std::make_pair(0.0f, 0.0f));  // default
	if (mit->size)
	  cptr->setDimensions(mit->size->width, mit->size->height);
	if (mit->symmetry)
	  cptr->setSymmetry(*(mit->symmetry));
	if (mit->site)
	  cptr->setSite(*(mit->site));

	lib.addCell(cptr);
      }
      if (cell_was_redefined && CheckPolicy<LEFERR_MACRO_REDEF>::abort)
	return false;
    }
    
    return true;
  }
  };
}
#endif

