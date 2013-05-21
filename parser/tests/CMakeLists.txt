foreach( tname basic lefbasic semantic )
  add_executable( ${tname} ${tname}.cpp )
  add_test( ${tname} ${tname} )
  if (MSVC)
    # fix compiler warning
    set_target_properties( ${tname} PROPERTIES COMPILE_FLAGS "/EHa" )
  endif()
  target_link_libraries( ${tname} lefdef )
endforeach()

add_executable( spefbasic spefbasic.cpp )
add_test( spefbasic spefbasic )
