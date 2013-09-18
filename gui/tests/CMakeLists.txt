# Qt5 support from KDAB requires component-level access at this time:
find_package( Qt5Widgets REQUIRED )

# automoc does not seem to work here for some reason, so:
QT5_WRAP_CPP(DS_MOC_SRC ../designscene.h)

add_executable ( basicgui basic.cpp ${DS_MOC_SRC} )
qt5_use_modules( basicgui Widgets )

if (MSVC)
  set_target_properties( basicgui PROPERTIES COMPILE_FLAGS "/EHa" )
endif()
add_test( BasicGUI basicgui )
if (WIN32)
  get_target_property( WidgetLibDir Qt5::Widgets LOCATION )
  get_filename_component( QT_BINARY_DIR ${WidgetLibDir} PATH )
  set_tests_properties( BasicGUI PROPERTIES ENVIRONMENT "Path=${QT_BINARY_DIR}" )
endif()



