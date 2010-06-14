// Tests for a Qt-based VLSI design viewer, part of EDASkel, a sample EDA app
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

// This file contains some simple unit tests for the GUI
#define BOOST_TEST_MODULE basic GUI tests
#include <boost/test/included/unit_test.hpp>

#include "../designscene.h"
#include "../../db/simpledb.h"

// Qt-using code requires a QApplication or QCoreApplication (if no GUI)
// apparently using QGraphicsScene requires QApplication even if no QGraphicsView
#include <QApplication>
#include <QGraphicsItem>

BOOST_AUTO_TEST_CASE( basic ) {
  int argc = 0;
  QApplication app(argc, 0);

  using namespace SimpleDB;
  // create a fake library
  Library lib;
  Library::CellPtr c1(new Library::Cell("BIGCELL"));
  c1->setDimensions(100, 100);
  lib.addCell(c1);
  Library::CellPtr c2(new Library::Cell("SMALLCELL"));
  c2->setDimensions(10, 10);
  lib.addCell(c2);

  // instantiate the library cells
  Database db;
  db.setDbuPerMicron(1);
  db.setExtent(Database::Rect(Database::Point(0, 0), Database::Point(10000, 10000)));
  Database::InstIter beg, end;
  tie(beg, end) = db.getInstances();
  BOOST_CHECK( beg == end );   // we haven't added any yet
  Database::InstPtr i1(new Instance("Inst1", "BIGCELL"));
  i1->setPlacement(Database::Point(2000, 2000), "N");
  db.addInst(i1);
  Database::InstPtr i2(new Instance("Inst2", "SMALLCELL"));
  i2->setPlacement(Database::Point(2095, 2095), "N");   // create a small overlap
  db.addInst(i2);
  tie(beg, end) = db.getInstances();
  BOOST_CHECK( distance(beg, end) == 2 );

  // now verify that the resulting scene has the expected contents
  DesignScene<Database, Library> myScene(db, lib);
  QRectF sr = myScene.sceneRect();
  // scene rectangle.  Remember that Y axis is reversed in Qt (0,0 is upper left)
  BOOST_CHECK( (sr.left() == 0.0) && (sr.right() == 10000.0) &&
	       (sr.bottom() == 10000.0) && (sr.top() == 0.0) );

  // check inserted items
  BOOST_CHECK( myScene.items().size() == 2 );
  QList<QGraphicsItem*> ilist = myScene.items();
  bool foundsmall = false, foundbig = false;
  while (!ilist.isEmpty()) {
    QGraphicsItem* i = ilist.takeFirst();
    QRectF br = i->boundingRect();
    if ((br.left() == 2000.0) && (br.right() == 2100.0) &&
	(br.bottom() == 2100.0) && (br.top() == 2000.0))
      foundbig = true;
    else if ((br.left() == 2095.0) && (br.right() == 2105.0) &&
	(br.bottom() == 2105.0) && (br.top() == 2095.0))
      foundsmall = true;
  }
  BOOST_CHECK( foundsmall && foundbig );

}
