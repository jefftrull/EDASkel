# Code in this directory requires Boost and in most cases also Eigen
find_package( Eigen3 NO_MODULE )
if( EIGEN3_FOUND )
  include_directories( ${CMAKE_SOURCE_DIR} )

  add_subdirectory( parasitics )
  add_subdirectory( tests )
else()
  message( "Eigen3 does not seem to be installed.  Not all executables will be built" )
endif()
