# Tcl Integration code.  Generic stuff needed for apps integrating with Tcl

# libraries: Qt
find_package(Qt4 4.6.2 REQUIRED)
include(${QT_USE_FILE})

# Tcl too
# BOZO figure out how to require 8.4 Tcl without messages about Tk
find_package(TCL 8.4)
include_directories(${TCL_INCLUDE_PATH})

QT4_WRAP_CPP( QTN_MOC_SRC qttclnotifier.h )
add_library(qttclnotifier qttclnotifier.cpp ${QTN_MOC_SRC})

