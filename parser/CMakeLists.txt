# if you want to debug Qi (the parser):
#add_definitions("-DBOOST_SPIRIT_DEBUG")
# and for the lexer
#add_definitions("-DBOOST_SPIRIT_LEXERTL_DEBUG")

# Phoenix v3 resolves some issues, and we should move to it anyway
add_definitions("-DBOOST_SPIRIT_USE_PHOENIX_V3")
# An issue with keywords, which should be fixed in 1.56, requires this:
add_definitions("-DBOOST_NO_CXX11_DECLTYPE_N3276")

add_subdirectory( tests )
add_library( lefdef lefdef.cpp )
