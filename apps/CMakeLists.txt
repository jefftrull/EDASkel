find_package(Qt4 4.6.2 REQUIRED)
include(${QT_USE_FILE})
# here's where we handle the Qt "moc" process
# have to make it a bit more explicit than in a qmake-based flow
QT4_WRAP_CPP(DS_MOC_SRC ../gui/designscene.h)
QT4_WRAP_CPP(DV_MOC_SRC ../gui/designview.h)
add_executable( sv simple_viewer.cpp ${DS_MOC_SRC} ${DV_MOC_SRC})
target_link_libraries( sv lefdef ${QT_LIBRARIES} ${Boost_LIBRARIES} )

