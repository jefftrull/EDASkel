cmake_minimum_required( VERSION 2.6 )
PROJECT( new_app )
ENABLE_TESTING()
find_package( Boost 1.43.0 REQUIRED )

add_subdirectory( parser )
add_subdirectory( gui )
add_subdirectory( db )
