// Semantic analysis for a Boost Spirit-based DEF parser, part of EDASkel, a sample EDA app
// Copyright (C) 2010 Jeffrey Elliot Trull <linmodemstudent@gmail.com>
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

#include "deftypes.h"

#if !defined(EDASKEL_DEF_SEMANTIC)
#define EDASKEL_DEF_SEMANTIC

namespace DefParse {

// we need to perform various checks, and then add DEF items to the database
// It's my feeling that we want to separate (most of) the checking from the parsing process
// This is less efficient, because in theory we should be able to blast design input directly
// into a database, but it simplifies and clarifies error handling.

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
		       DEFERR_SITE_OUTSIDE_DIEAREA,
		       DEFERR_DIE_MALFORMED_RECT,
		       DEFERR_INST_UNKNOWN_CELL,
		       DEFERR_INST_OVERLAP,
		       DEFERR_INST_OUTSIDE_DIEAREA
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
  // Look at simpledb.h for examples of the "Concept" DB must implement
  template<class DB, class Lib, template<DefCheckError> class CheckPolicy = DefCheckPolicy>
    class DefChecker {
  public:
    bool CheckAndInsert(const def& defin, const Lib& lib, DB& db) {
      typedef typename DB::Rect Rect;
      typedef typename DB::Point Point;

      // diearea becomes design extent
      // check for diearea is zero or negative
      bool die_malformed_rect = ((defin.diearea.ll.x >= defin.diearea.ur.x) ||
				 (defin.diearea.ll.y >= defin.diearea.ur.y));
      if (die_malformed_rect && !CheckPolicy<DEFERR_DIE_MALFORMED_RECT>::silent) {
	std::cerr << "DEF Checker: malformed die area (" << defin.diearea.ll.x << " " << defin.diearea.ll.y;
	std::cerr << ") (" << defin.diearea.ur.x << " " << defin.diearea.ur.y << ")\n";
      }
      if (!die_malformed_rect || !CheckPolicy<DEFERR_DIE_MALFORMED_RECT>::skip) {
	db.setExtent(Rect(Point(defin.diearea.ll.x, defin.diearea.ll.y),
			  Point(defin.diearea.ur.x, defin.diearea.ur.y)));
      }
      if (die_malformed_rect && CheckPolicy<DEFERR_DIE_MALFORMED_RECT>::abort) {
	return false;
      }

      // Examine site statements and see if they make sense (exist and are inside design)

      // Examine each instance
    
      return true;
    }
  };
}
#endif

