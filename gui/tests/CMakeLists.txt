QT4_WRAP_CPP(DS_MOC_SRC ../designscene.h)
add_executable ( basicgui basic.cpp ${DS_MOC_SRC})
add_test( BasicGUI basicgui )
target_link_libraries(basicgui ${QT_LIBRARIES})


