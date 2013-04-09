find_package( Qt5Widgets REQUIRED )

# here's where we handle the Qt "moc" process
# have to make it a bit more explicit than in a qmake-based flow
QT5_WRAP_CPP(DS_MOC_SRC ../gui/designscene.h)
QT5_WRAP_CPP(DV_MOC_SRC ../gui/designview.h)
add_executable( sv simple_viewer.cpp ${DS_MOC_SRC} ${DV_MOC_SRC})
qt5_use_modules( sv Widgets )

target_link_libraries( sv lefdef ${Boost_LIBRARIES} )
