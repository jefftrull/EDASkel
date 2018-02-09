find_package( Qt5Widgets REQUIRED )

# here's where we handle the Qt "moc" process
# have to make it a bit more explicit than in a qmake-based flow
QT5_WRAP_CPP(DV_MOC_SRC ../gui/designview.h)

if( MSVC )
  # prevents linker errors (arg, default_line_length)
  add_definitions( -DBOOST_PROGRAM_OPTIONS_DYN_LINK )
endif()

add_executable( spef load_spef.cpp )
target_link_libraries( spef Boost::boost )

include_directories(${TCL_INCLUDE_PATH})

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  # suppress issue in Qt5 headers
  add_definitions( -Wno-inconsistent-missing-override )
endif()

if( TCL_FOUND )
  add_executable( sv simple_viewer.cpp ${DV_MOC_SRC})
  target_link_libraries( sv Qt5::Widgets lefdef qttclnotifier ${TCL_LIBRARY} Boost::program_options )
endif()

