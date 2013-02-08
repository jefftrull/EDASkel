cmake_minimum_required( VERSION 2.8 )
PROJECT( new_app )
ENABLE_TESTING()
find_package( Boost 1.52 REQUIRED COMPONENTS program_options )
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

# global flags
if (CMAKE_COMPILER_IS_GNUCC)
  add_definitions(-Wall)          # be conservative about checks
  add_definitions(-std=c++0x)     # enable C++11
endif()

add_subdirectory( parser )
add_subdirectory( gui )
add_subdirectory( db )
add_subdirectory( apps )
add_subdirectory( tclint )
