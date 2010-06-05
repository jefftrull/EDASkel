// A remedial design database and cell library, part of EDASkel, a sample EDA app
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

#ifndef SIMPLEDB_H
#define SIMPLEDB_H

// the idea here is to provide a simple database that conforms to the "Concept" expected
// by other parts of EDASkel.  It won't be fast, but should be easy to understand.
// Wrappers for "real" design databases are preferable for production work.

#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include <boost/assert.hpp>
#include <boost/optional.hpp>
// using the TR1, not Boost, versions of things wherever possible
// (I like Boost, but standards are even better)
#include <tr1/memory>

namespace SimpleDB {

  class LibCell {
  public:
    LibCell(const std::string& name) : m_name(name) {}
    std::string getName() const { return m_name; }

    void setDimensions(int w, int h) {m_dimensions = std::make_pair(w, h);}
    bool hasDimensions() const {return m_dimensions;}
    int getWidth() const {BOOST_ASSERT(m_dimensions); return m_dimensions->first;}
    int getHeight() const {BOOST_ASSERT(m_dimensions); return m_dimensions->second;}

  private:
    std::string m_name;
    // I like to use boost::optional because it packages the idea of "data not present" cleanly
    // We may have taken this cell from a library that doesn't have physical dimensions for the cell
    // (e.g. a timing or power library) or it may not have been specified for some reason
    boost::optional<std::pair<int, int> > m_dimensions;
  };

  class Library {
  public:
    typedef LibCell Cell;
    typedef std::tr1::shared_ptr<Cell> CellPtr;
    Library() {}

    CellPtr findCell(const std::string& name) const {
      // linear search for a matching name
      // could make this Phoenix or Lambda or a complex bind expression.  But should I?
      for (std::vector<CellPtr>::const_iterator cit = m_cells.begin();
	   cit != m_cells.end(); ++cit) {
	if ((*cit)->getName() == name)
	  return *cit;
      }
      return CellPtr(); // the shared_ptr equivalent of NULL
    }
					     
    void addCell(CellPtr cell) {
      m_cells.push_back(cell);
    }

  private:
    std::vector<CellPtr> m_cells;
  };

  class DesPoint {
  public:
  DesPoint(int x, int y) : m_x(x), m_y(y) {}
    int x() const {return m_x;}
    int y() const {return m_y;}
  private:
    int m_x, m_y;
  };
    

  class DesRect {
  public:
  DesRect(DesPoint ll, DesPoint ur) : m_ll(ll), m_ur(ur) {}
    DesPoint ll() const {return m_ll;}
    DesPoint ur() const {return m_ur;}
  private:
    DesPoint m_ll, m_ur;
  };

  class Instance {
  public:
  Instance(const std::string& name, const std::string& cname) : m_name(name), m_cname(cname) {}
    // ultimately this will be a pointer to a library cell
    std::string getCellName() const {return m_cname;}

    std::string getName() const {return m_name;}

    void setOrigin(int x, int y) {m_origin = std::make_pair(x, y);}
    bool hasOrigin() const {return m_origin;}
    DesPoint getOrigin() const {BOOST_ASSERT(m_origin);
                                return DesPoint(m_origin->first, m_origin->second);}
  private:
    std::string m_name, m_cname;
    boost::optional<std::pair<int, int> > m_origin;
  };

  class Database {
  public:
    typedef DesPoint Point;
    typedef DesRect Rect;
    typedef std::tr1::shared_ptr<Instance> InstPtr;
    typedef std::vector<InstPtr>::const_iterator InstIter;  // users cannot change stored pointers

    Database() {};
    std::pair<InstIter, InstIter> getInstances() const { return std::make_pair(m_instances.begin(),
									       m_instances.end()); }
    void setExtent(DesRect r) {m_extent = r;}
    bool hasExtent() const {return m_extent;}
    DesRect getExtent() const {BOOST_ASSERT(m_extent); return *m_extent;}

    void addInst(InstPtr i) {
      m_instances.push_back(i);
    }

  private:
    std::vector<InstPtr> m_instances;
    boost::optional<DesRect> m_extent;
  };

}    

#endif
