# Qt5 support from KDAB requires component-level access at this time:
find_package( Qt5Widgets REQUIRED )

add_executable ( basicgui basic.cpp )
target_link_libraries( basicgui Qt5::Widgets Boost::boost )

if (MSVC)
  set_target_properties( basicgui PROPERTIES COMPILE_FLAGS "/EHa" )
endif()
add_test( BasicGUI basicgui )
if (WIN32)
  get_target_property( WidgetLibDir Qt5::Widgets LOCATION )
  get_filename_component( QT_BINARY_DIR ${WidgetLibDir} PATH )
  set_tests_properties( BasicGUI PROPERTIES ENVIRONMENT "Path=${QT_BINARY_DIR}" )
endif()



