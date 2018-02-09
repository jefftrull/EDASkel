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
#include <QGraphicsItem>
#include <QPainter>

namespace EDASkel {

class DesignScene : public QGraphicsScene {
  Q_OBJECT

  // base class of templated child class so we can "erase" the db/lib details and use
  // the highlight interface uniformly
  class InstItemBase : public QGraphicsItem {
    bool highlighted_;
  public:
    InstItemBase() : highlighted_(false)  {}
    ~InstItemBase() = default;

    void setHighlight(bool highlighted) {
      highlighted_ = highlighted;
      update(boundingRect());
    }

    virtual QRectF boundingRect() const = 0;

    void paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*) {
      QPen pen = highlighted_ ? QPen(QColor::fromRgb(255, 215, 0), 0) // golden
                              : QPen(Qt::red, 0);
      painter->setPen(pen);
      QRectF boundary = boundingRect();
      painter->drawRect(boundary);
      // add dogear to indicate origin corner (lower left, for N instances)
      painter->drawLine(QPointF(0, boundary.bottom()/2),
                        QPointF(boundary.right()/2, 0));
    }

  };

  template<class DB, class Lib>
  class InstItem : public InstItemBase {
    typename DB::InstPtr inst_;
    typename Lib::CellPtr cell_;

  public:
    InstItem(typename DB::InstPtr inst, typename Lib::CellPtr cell) : inst_(inst), cell_(cell) {
      // handle orientation
      // It looks like W, S, E = 90, 180, -90 degree rotations
      // The "flip" versions can happen by mirroring around the Y axis first
      std::string ort = inst_->getOrient();
      // right-hand rotation (natural, positive angle) for each orientation
      std::map<std::string, int> rotations{{"N", 0}, {"W", 90}, {"S", 180}, {"E", -90}};

      QTransform xform;
      if (ort.front() == 'F') {
        ort = ort.back();
        xform.scale(-1.0, 1.0);               // apply flip
        xform.rotate(-rotations.at(ort));     // negate rotation to compensate for flip
      } else {
        // Qt's clockwise rotation is cancelled out by the fact that the Z axis is into
        // the screen due to our overall flipping in the "view".  So use the natural value.
        xform.rotate(rotations.at(ort));
      }

      // Finally, the position specified in DEF is that of the lower left *after* orientation
      // figure out where the new LL is so we can compensate
      QPointF newLL = xform.mapRect(boundingRect()).topLeft();  // "top" left due to scene flip
      xform = xform * QTransform::fromTranslate(-newLL.x(), -newLL.y());

      setTransform(xform);

      // create a formatted "Tool Tip" (hover text) to identify this inst
      typename DB::Point orig = inst_->getOrigin();
      setToolTip(QStringLiteral("%1 (%2) (%3, %4) %5").
                 arg(inst_->getName().c_str()).
                 arg(inst_->getCellName().c_str()).
                 arg(orig.x()).
                 arg(orig.y()).
                 arg(inst_->getOrient().c_str()));
    }
    
    QRectF boundingRect() const override {
      // look up cell dimensions
      // TODO: handle rectilinear polygon boundaries
      // TODO: handle other coordinate sizes (e.g., long, float)
      float width = cell_->getWidth();
      float height = cell_->getHeight();

      return QRectF(0, 0, width, height);
    }

  };

  // raw pointers OK here because this class owns (in Qt parlance is parent of) the InstItems
  InstItemBase* highlightedInst;
  std::map<std::string, InstItemBase*> instNameMap;

 public:
  ~DesignScene() = default;
  template<class DB, class Lib>
  explicit DesignScene(const DB& db,
		       const Lib& lib,
		       QObject* parent = nullptr)
      : QGraphicsScene(parent), highlightedInst(nullptr) {
    // basic setup
    setBackgroundBrush(Qt::black);
    int dbu = db.getDbuPerMicron();
    // assume for now that LEF is in microns

    // get design extent from DB and set scene accordingly
    typename DB::Rect design_extent = db.getExtent();
    // Qt uses "topleft/bottomright" which conveniently is exactly what I consider bottomleft/topright
    setSceneRect(QRectF(QPointF(design_extent.ll().x(), design_extent.ll().y()),
			QPointF(design_extent.ur().x(), design_extent.ur().y())));

    // get instances and add their proxy objects
    for (auto const& instp : db.getInstances()) {
      const typename Lib::CellPtr cell = lib.findCell(instp->getCellName());
      if (!cell)
	continue;   // or produce an error somehow?

      InstItemBase* inst = new InstItem<DB, Lib>(instp, cell);
      // InstItem works in LEF units so scale/apply position here
      inst->setTransform(inst->transform() * QTransform::fromScale(dbu, dbu));
      inst->setPos(instp->getOrigin().x(), instp->getOrigin().y());
      addItem(inst);
      instNameMap.emplace(instp->getName(), inst);
    }
  }

  bool highlightInstance(std::string const& name) {
    if (highlightedInst != nullptr) {
      highlightedInst->setHighlight(false);
    }
    auto map_it = instNameMap.find(name);
    if (map_it != instNameMap.end()) {
      highlightedInst = map_it->second;
      highlightedInst->setHighlight(true);
      return true;
    }
    highlightedInst = nullptr;
    return false;
  }

};

}

#endif
