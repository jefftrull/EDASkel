# if you want to debug Qi (the parser):
#add_definitions("-DBOOST_SPIRIT_DEBUG")
# and for the lexer
#add_definitions("-DBOOST_SPIRIT_LEXERTL_DEBUG")

add_subdirectory( tests )
if( MSVC )
  # in debug mode we can get C1128; this is the workaround
  set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj" )
endif()

add_library( lefdef lefdef.cpp )
