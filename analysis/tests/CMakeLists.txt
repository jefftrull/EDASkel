foreach( tname mor )
  add_executable( ${tname} ${tname}.cpp )
  add_test( ${tname} ${tname} )
endforeach()
