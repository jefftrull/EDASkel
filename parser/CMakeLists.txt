# if you want to debug Qi (the parser):
#add_definitions("-DBOOST_SPIRIT_DEBUG")
# and for the lexer
add_definitions("-DBOOST_SPIRIT_LEXERTL_DEBUG")

add_subdirectory( tests )
add_library( lefdef lefdef.cpp )
