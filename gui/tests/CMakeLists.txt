QT4_WRAP_CPP(DS_MOC_SRC ../designscene.h)
add_executable ( basicgui basic.cpp ${DS_MOC_SRC})
if (MSVC)
  set_target_properties( basicgui PROPERTIES COMPILE_FLAGS "/EHa" )
endif()
add_test( BasicGUI basicgui )
target_link_libraries(basicgui ${QT_LIBRARIES})


