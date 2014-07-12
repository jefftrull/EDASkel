foreach( fname coupling_eigen_odeint rlc_eigen_odeint
               coupling_odeint rlc_odeint coupling_prima_eigen_odeint coupling_sprim_eigen_odeint )
  add_executable( ${fname} ${fname}.cpp )
endforeach()

