// Semantic analysis for a Boost Spirit-based DEF parser, part of EDASkel, a sample EDA app
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

#include <boost/numeric/conversion/cast.hpp>

#include "deftypes.h"

#if !defined(EDASKEL_DEF_SEMANTIC)
#define EDASKEL_DEF_SEMANTIC

namespace DefParse {

// we need to perform various checks, and then add DEF items to the database
// It's my feeling that we want to separate (most of) the checking from the parsing process
// This is less efficient, because in theory we should be able to blast design input directly
// into a database, but it simplifies and clarifies error handling.  It also allows
// us to handle default values more cleanly

// Error Policies
// LEF/DEF is a funny sort of "standard".  Many existing tools accept input that is theoretically
// wrong, and many user flows have grown to expect that their idiosyncratic input will work.
// Therefore we need to provide some flexibility in how the parser handles errors.

// I want the user to be able to specify for each kind of error, separately whether to
// produce an error or remain silent, insert anyway (where possible) or not, and
// abort the entire semantic check or not.  I also want to supply reasonable defaults,
// and have it all handled by the compiler (i.e., not check these policies at runtime).

// I think the cleanest approach is to specify a single policy for error handling that users
// can override.  If we use an enum to represent different kinds of errors, I can give the generic
// implementation for all values and users can "specialize" for whichever errors they want

  // list of known errors
  enum DefCheckError { DEFERR_SITE_UNKNOWN_NAME,
		       DEFERR_SITE_MALFORMED,
		       DEFERR_SITE_OUTSIDE_DIEAREA,
		       DEFERR_DIEAREA_MALFORMED,
		       DEFERR_INST_UNKNOWN_CELL,
		       DEFERR_INST_REDEFINED
  };

  // default handling approach
  template<DefCheckError>
  class DefCheckPolicy {
  public:
    static bool const silent = false;  // do not be silent (give an error message)
    static bool const skip = true;     // do skip (i.e., do not take action on this element)
    static bool const abort = false;   // do not abort (i.e., try to keep going)
    // note: skip is considered before abort, so you can insert something and then abort anyway
  };

  // and now a member function to check the results of a DEF parse and conditionally insert
  // the elements into a database
  // Look at simpledb.hpp for examples of the "Concept" DB must implement
  template<class DB, class Lib, template<DefCheckError> class CheckPolicy = DefCheckPolicy>
    class DefChecker {
  public:
    bool CheckAndInsert(const def& defin, const Lib& lib, DB& db) {
      using boost::numeric_cast;

      typedef typename DB::Rect Rect;
      typedef typename DB::Point Point;

      // diearea becomes design extent
      // check for diearea is zero or negative
      bool die_malformed_rect = ((defin.diearea.ll.x >= defin.diearea.ur.x) ||
				 (defin.diearea.ll.y >= defin.diearea.ur.y));
      if (die_malformed_rect && !CheckPolicy<DEFERR_DIEAREA_MALFORMED>::silent) {
	std::cerr << "DEF Checker: malformed die area (" << defin.diearea.ll.x << " " << defin.diearea.ll.y;
	std::cerr << ") (" << defin.diearea.ur.x << " " << defin.diearea.ur.y << ")\n";
      }
      if (!die_malformed_rect || !CheckPolicy<DEFERR_DIEAREA_MALFORMED>::skip) {
	db.setExtent(Rect(Point(defin.diearea.ll.x, defin.diearea.ll.y),
			  Point(defin.diearea.ur.x, defin.diearea.ur.y)));
      }
      if (die_malformed_rect && CheckPolicy<DEFERR_DIEAREA_MALFORMED>::abort) {
	return false;
      }

      // verify units were supplied?
      // also: only certain combinations of dbupermicron in DEF vs. LEF are allowed - check this
      db.setDbuPerMicron(defin.dbupermicron);

      // Examine row/site statements and see if they make sense (exist and are inside design)
      // if STEP is provided at least one must be 1
      for (auto const& rs : defin.rows) {
	// verify the Site referenced by this DEF line is defined in the library
	// (or tech LEF, tbd)
	typedef typename Lib::SitePtr LibSitePtr;
	LibSitePtr sptr = lib.findSite(rs.sitename);
	bool site_not_found = (sptr == LibSitePtr());
	if (site_not_found && !CheckPolicy<DEFERR_SITE_UNKNOWN_NAME>::silent)
	  std::cerr << "DEF Checker: unknown site name " << rs.sitename << std::endl;

	// set up some defaults for count and stepping
	int xcount = 1; int ycount = 1;
	int xstep = 0; int ystep = 0;
	if (!site_not_found) {
	  // "dbupermicron" is expected to be large enough so these calculations produce an integer result:
          xstep = numeric_cast<int>(numeric_cast<float>(defin.dbupermicron) * sptr->getWidth());
          ystep = numeric_cast<int>(numeric_cast<float>(defin.dbupermicron) * sptr->getHeight());
	}
	if (rs.repeat) {
	  xcount = rs.repeat->xrepeat;
	  ycount = rs.repeat->yrepeat;
	  if (rs.repeat->step) {
	    xstep = rs.repeat->step->first;
	    ystep = rs.repeat->step->second;
	  }
	}
	bool site_malformed = (ycount != 1) && (xcount != 1);
	if (site_malformed && !CheckPolicy<DEFERR_SITE_MALFORMED>::silent) {
	  if (rs.rowname)
	    std::cerr << "DEF Checker: row " << *rs.rowname << " must have either X or Y count == 1\n";
	  else
	    std::cerr << "DEF Checker: site row definition must have either X or Y count == 1\n";
	}

	// check limits
	int farx = rs.origx + (xcount - 1) * xstep;
	int fary = rs.origy + (ycount - 1) * ystep;
	if (!site_not_found) {
	  // also add in site width/height
          farx += numeric_cast<int>(numeric_cast<float>(defin.dbupermicron) * sptr->getWidth());
          fary += numeric_cast<int>(numeric_cast<float>(defin.dbupermicron) * sptr->getHeight());
	}
	bool out_of_diearea = ((rs.origx < defin.diearea.ll.x) || (farx < defin.diearea.ll.x) ||
			       (rs.origy < defin.diearea.ll.y) || (fary < defin.diearea.ll.y) ||
			       (rs.origx > defin.diearea.ur.x) || (farx > defin.diearea.ur.x) ||
			       (rs.origy > defin.diearea.ur.y) || (fary > defin.diearea.ur.y));
	if (out_of_diearea && !CheckPolicy<DEFERR_SITE_OUTSIDE_DIEAREA>::silent) {
	  if (rs.rowname)
	    std::cerr << "DEF Checker: row " << *rs.rowname << " from (" << rs.origx << ", " << rs.origy << ") to (" << farx << ", " << fary << ") extends beyond DIEAREA boundary\n";
	  else
	    std::cerr << "DEF Checker: sites from (" << rs.origx << ", " << rs.origy << ") to (" << farx << ", " << fary << ") extend beyond DIEAREA boundary\n";
	}
	  
	if ((site_not_found && CheckPolicy<DEFERR_SITE_UNKNOWN_NAME>::skip) ||
	    (out_of_diearea && CheckPolicy<DEFERR_SITE_OUTSIDE_DIEAREA>::skip) ||
	    (site_malformed && CheckPolicy<DEFERR_SITE_MALFORMED>::skip))
	  continue;
	// add sites to database
	if (rs.rowname)
	  db.addRow(*(rs.rowname), rs.sitename, xcount, ycount, Point(rs.origx, rs.origy), xstep, ystep);
	else
	  db.addRow("<UNNAMED>", rs.sitename, xcount, ycount, Point(rs.origx, rs.origy), xstep, ystep);
	if ((site_not_found && CheckPolicy<DEFERR_SITE_UNKNOWN_NAME>::abort) ||
	    (out_of_diearea && CheckPolicy<DEFERR_SITE_OUTSIDE_DIEAREA>::abort) ||
	    (site_malformed && CheckPolicy<DEFERR_SITE_MALFORMED>::abort))
	  return false;

      }

      // Examine each instance
      for (auto const& inst : defin.components) {
	typedef typename Lib::CellPtr CellPtr;
	CellPtr cell = lib.findCell(inst.celltype);
	bool cell_not_found = (cell == CellPtr());
	if (cell_not_found && !CheckPolicy<DEFERR_INST_UNKNOWN_CELL>::silent)
	  std::cerr << "DEF Checker: Instance " << inst.name << " is of unknown cell type " << inst.celltype << std::endl;

	typedef typename DB::InstPtr InstPtr;
	InstPtr previnst = db.findInst(inst.name);
	bool inst_redefined = (previnst != InstPtr());
	if (inst_redefined && !CheckPolicy<DEFERR_INST_REDEFINED>::silent)
	  std::cerr << "DEF Checker: Instance " << inst.name << " is of unknown cell type " << inst.celltype << std::endl;

	if ((!cell_not_found || !CheckPolicy<DEFERR_INST_UNKNOWN_CELL>::skip) &&
	    (!inst_redefined || !CheckPolicy<DEFERR_INST_REDEFINED>::skip)) {
	  if (inst_redefined)
	    // we really can't have both of these in the db, so delete the previously found one
	    db.removeInst(previnst);
	  typedef typename DB::Inst Inst;
	  InstPtr instp(new Inst(inst.name, inst.celltype));
	  if (inst.placement)
	    instp->setPlacement(Point(inst.placement->origin.x, inst.placement->origin.y),
			       inst.placement->orient, (inst.placement->plcfix == "FIXED"));
	  db.addInst(instp);
	}
	if ((cell_not_found && CheckPolicy<DEFERR_INST_UNKNOWN_CELL>::abort) ||
	    (inst_redefined && CheckPolicy<DEFERR_INST_REDEFINED>::abort))
	  return false;
      }
    
      return true;
    }
  };
}
#endif

