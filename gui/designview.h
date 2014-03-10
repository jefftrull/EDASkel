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

/* -*-C++-*- */

// This is the "view" part of the Qt Model/View pattern, applied to designs

#ifndef EDASKEL_DESIGNVIEW_H
#define EDASKEL_DESIGNVIEW_H

#include <memory>

#include <QGraphicsView>
#include <QKeyEvent>

namespace EDASkel {

class DesignView : public QGraphicsView {
  Q_OBJECT
public:
  explicit DesignView(std::unique_ptr<QGraphicsScene> scene,
                      QWidget* parent = 0 ) : QGraphicsView(scene.get(), parent),
                                              scene_(std::move(scene)) {

    // compensate for Qt's use of 0,0 to mean "upper left"
    // as designers we expect to see coordinates increasing up and to the right
    scale(1.0, -1.0);  // the scene is upside down in Qt coordinates

  }
  void keyPressEvent(QKeyEvent* event) {
    // intercept a few keystrokes so we can zoom and "fit"
    if (event->key() == Qt::Key_F)
      fitInView(sceneRect(), Qt::KeepAspectRatio);  // show all of design
    else if (event->key() == Qt::Key_Q)
      close();                                      // close viewer (and quit)
    else if (event->key() == Qt::Key_I)
      scale(1.5, 1.5);                              // zoom in
    else if (event->key() == Qt::Key_O)
      scale(0.67, 0.67);                            // zoom out
    else
      QGraphicsView::keyPressEvent(event);          // just pass it on
  }

private:
  std::unique_ptr<QGraphicsScene> scene_;

};

}
#endif
