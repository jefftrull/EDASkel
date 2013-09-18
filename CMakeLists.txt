cmake_minimum_required( VERSION 2.8.9 )

PROJECT( new_app )
ENABLE_TESTING()
find_package( Boost 1.53 REQUIRED COMPONENTS program_options )
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
if (MSVC)
  # prevent auto linking
  add_definitions( -DBOOST_ALL_NO_LIB=1 )
endif()

# global flags
if (CMAKE_COMPILER_IS_GNUCC)
  add_definitions(-Wall)          # be conservative about checks
  add_definitions(-std=c++0x)     # enable C++11
elseif (MSVC11)
  # VC11 workaround: emplace_back requires variadic template emulation, for which VC11 has
  # a default limit of 5 arguments due to some memory explosion issue.  Raising this to 6
  # allows our calls to emplace_back (to name one example) to succeed
  add_definitions(-D_VARIADIC_MAX=6)
endif()

# for packages without a standard CMake find script:
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

add_subdirectory( parser )
add_subdirectory( gui )
add_subdirectory( db )
add_subdirectory( apps )
add_subdirectory( tclint )
add_subdirectory( analysis )
