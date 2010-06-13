include_directories( ${Boost_INCLUDE_DIRS} )
add_subdirectory( tests )
add_library( lefdef lefdef.cpp )
