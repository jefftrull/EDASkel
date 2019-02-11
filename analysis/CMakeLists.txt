# Code in this directory requires Boost and in most cases also Eigen
find_package( Eigen3 )
if( EIGEN3_FOUND )
  include_directories( ${CMAKE_SOURCE_DIR} )
  if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Clang warns on some Eigen things so we need to do this:
    add_definitions( -Wno-sign-compare -Wno-deprecated-register -Wno-deprecated-declarations )
  endif()

  add_subdirectory( parasitics )
  add_subdirectory( tests )
else()
  message( "Eigen3 does not seem to be installed.  Not all executables will be built" )
endif()
