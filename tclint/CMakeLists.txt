# Tcl Integration code.  Generic stuff needed for apps integrating with Tcl

# libraries: Qt (Core only)
find_package( Qt5Core )

include_directories(${TCL_INCLUDE_PATH})

QT5_WRAP_CPP( QTN_MOC_SRC qttclnotifier.h )
add_library(qttclnotifier qttclnotifier.cpp ${QTN_MOC_SRC})
qt5_use_modules( qttclnotifier Core )
