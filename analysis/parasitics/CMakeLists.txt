# Code in this directory requires Boost and in most cases also Eigen
foreach( fname coupling_odeint rlc_odeint )
  add_executable( ${fname} ${fname}.cpp )
endforeach()

find_package( Eigen3 )
if( EIGEN3_FOUND )
  include_directories( ${EIGEN3_INCLUDE_DIR} )
  foreach( fname coupling_eigen_odeint rlc_eigen_odeint coupling_prima_eigen_odeint
                 rlc_prima_eigen_odeint )
    add_executable( ${fname} ${fname}.cpp )
  endforeach()
else()
  message( "Eigen3 does not seem to be installed.  Not all executables will be built" )
endif()

