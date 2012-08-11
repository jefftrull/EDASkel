// A Qt-based VLSI design viewer, part of EDASkel, a sample EDA app
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

// Qt provides a model/view based way of displaying user data where you can use
// the user coordinate system and describe objects in a natural way (the model), and
// separately use Qt's QGraphicsView classes to "view" your data.  These built-in
// view classes take the burden of handling most user interaction (zoom, pan, etc.)
// away from the programmer and allow you to concentrate on what your data is doing.

// this class implements a Qt "Scene" (i.e., a model) based on templated database and
// library classes.  Programmers who want to use a particular database may have to write
// a thin wrapper class that implements the required methods and typedefs.  Doing it this
// way provides agnosticism w.r.t the database, while providing some performance improvement
// vs. C-style callbacks via function pointers, because the target database API will be
// visible to the compiler.

/* -*-C++-*- */
#ifndef EDASKEL_DESIGNSCENE_H
#define EDASKEL_DESIGNSCENE_H


#include <QGraphicsScene>
#include <QGraphicsRectItem>

namespace EDASkel {

// Templated classes can't be turned directly into QObjects by "moc" so I have to
// separate DB/Lib access stuff into a child class and leave Signal/Slot stuff
// in a non-templated base class.
class DesignSceneBase : public QGraphicsScene {
  Q_OBJECT
  public:
  explicit DesignSceneBase( QObject* parent = 0 ) : QGraphicsScene(parent) {}
  ~DesignSceneBase() {}
};

template<class DB, class Lib>
class DesignScene : public DesignSceneBase {
 public:
  explicit DesignScene(const DB& db,
		       const Lib& lib,
		       QObject* parent = 0) : DesignSceneBase(parent) {
    // basic setup
    setBackgroundBrush(Qt::black);
    int dbu = db.getDbuPerMicron();
    // assume for now that LEF is in microns

    // get design extent from DB and set scene accordingly
    typename DB::Rect design_extent = db.getExtent();
    // Qt uses "topleft/bottomright" which conveniently is exactly what I consider bottomleft/topright
    setSceneRect(QRectF(QPointF(design_extent.ll().x(), design_extent.ll().y()),
			QPointF(design_extent.ur().x(), design_extent.ur().y())));
    // get instances and add their boundaries
    // ultimately you want to subclass one of the GraphicsItem classes instead
    // so you can store a pointer to the instance and maybe further customize
    // the appearance (e.g., by marking the UR corner to make orientation more obvious)
    QPen instPen(Qt::red);
    typename DB::InstIter iit, end;
    for (std::tie(iit, end) = db.getInstances(); iit != end; ++iit) {
      const typename Lib::CellPtr cell = lib.findCell((*iit)->getCellName());
      if (!cell)
	continue;   // or produce an error somehow?

      // look up cell dimensions
      // TODO: handle rectilinear polygon boundaries
      // TODO: handle other coordinate sizes (e.g., long, float)
      int width = cell->getWidth() * dbu;
      int height = cell->getHeight() * dbu;

      typename DB::Point orig = (*iit)->getOrigin();
      // handle orientation
      // From what I can tell the DEF coordinate is the location of the LL corner *after*
      // orientation is taken into account
      // For instance outlines all we care about is the location of the UR corner
      // It seems like this can be handled by conditionally exchanging width and height
      if (((*iit)->getOrient() == "E") || ((*iit)->getOrient() == "FE") ||
	  ((*iit)->getOrient() == "W") || ((*iit)->getOrient() == "FW"))
	std::swap(width, height);
      // when we display the geometries within the cell we'll have to revisit this with
      // a more sophisticated technique involving transforms

      QGraphicsRectItem* instrect = addRect(QRectF(orig.x(), orig.y(), width, height));
      instrect->setPen(instPen);
      // create a formatted "Tool Tip" (hover text) to identify this inst
      instrect->setToolTip(QString("%1 (%2) (%3, %4) %5").
			   arg((*iit)->getName().c_str()).
			   arg((*iit)->getCellName().c_str()).
			   arg(orig.x()).
			   arg(orig.y()).
			   arg((*iit)->getOrient().c_str()));
    }
  }
};

}

#endif
