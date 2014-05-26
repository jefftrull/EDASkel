// Test of odeint using intersignal coupling
// Now with the use of Eigen to separate the equations
// and to perform PRIMA model reduction
// Jeff Trull 2013-06-04

#include <vector>
#include <numeric>
using namespace std;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric;
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;
 
#include "analysis/mna.hpp"
#include "analysis/prima.hpp"
#include <Eigen/SparseQR>

// Utility functions for working with Eigen
template<class M>
bool canLDLTDecompose(const M& m) {
   // m must be positive or negative semidefinite, which means its eigenvalues
   // must all be real-valued, and either all non-positive or all non-negative.
   // Use the magic of Eigen reductions to implement:
   auto eigenvalues = EigenSolver<M>(m).eigenvalues();
   return (eigenvalues.array().imag() == 0.0).all() &&   // all real
      ((eigenvalues.array().real() >= 0.0).all() ||      // non-negative
       (eigenvalues.array().real() <= 0.0).all());       // or non-positive
}

typedef vector<double> state_type;

struct signal_coupling {
  static const size_t q = 8; // desired number of state variables in the reduced system
                             // 12 is the "natural" count (1 per circuit node)
  static const size_t N = 3; // port count (one per input/output: two drivers and the victim rcvr)

  // final system for simulation.  State variable count determined by regularization:
  Matrix<double, Dynamic, Dynamic> coeff_;    // reduced system state evolution
  Matrix<double, Dynamic, 2> input_;          // inputs (agg/vic) to reduced system state
  Matrix<double, 1, Dynamic> output_;         // reduced system state to chosen outputs

  // original MNA matrices: 12 circuit nodes, 3 independent sources
  typedef Matrix<double, 15, Dynamic> Matrix15dX;
  Matrix15dX Xfinal;

  // waveform control
  double agg_slew_, agg_start_;
  double vic_slew_, vic_start_;
  double v_;

  signal_coupling(
    double agg_r1,             // aggressor first stage pi model (prior to coupling point)
    double agg_c1,
    double agg_r2,             // aggressor second stage pi model (after coupling point)
    double agg_c2,
    double agg_cl,             // aggressor final load cap
    double agg_slew,           // aggressor slew rate (V/s)
    double agg_imp,            // aggressor driver impedance (placed after voltage source)
    double agg_start,          // aggressor driver start time (the *center* of the ramp!)

    double vic_r1,             // victim first stage pi model (prior to coupling point)
    double vic_c1,
    double vic_r2,             // victim second stage pi model (after coupling point)
    double vic_c2,
    double vic_cl,             // victim final load cap
    double vic_slew,           // victim input slew rate (V/s)
    double vic_imp,            // victim driver impedance
    double vic_start,          // victim driver start time

    double coup_c,             // coupling capacitance placed at central point
    double v                   // power supply (and thus max) voltage
    ) :
    agg_slew_(agg_slew), agg_start_(agg_start),
    vic_slew_(vic_slew), vic_start_(vic_start),
    v_(v) {
      // Use Eigen to construct the matrix we will use to calculate the dV/dt values.

      // Formerly I wrote out KCL by hand for every node in the circuit, then solved
      // manually to put the equations in the form needed by odeint.  With a matrix
      // library such as Eigen, and the Modified Nodal Analysis technique, this becomes
      // much simpler.  Each resistor and capacitor has a "stamp" of conductance values
      // to sum into one of two matrices which multiply the state or state derivative
      // matrices, respectively.  This produces the matrix equation C*dV/dt = -G*V,
      // which is equivalent to the KCL equations.  This can be solved for dV/dt in terms
      // of V, giving us the desired format for odeint.

      // apply values via "stamp"
      std::vector<Eigen::Triplet<double> > Centries, Gentries;

      using namespace EDASkel::analysis::mna;
      stamp(Gentries, 0, 1, 1/agg_imp);   // driver impedances
      stamp(Gentries, 6, 7, 1/vic_imp);

      stamp(Centries, 1, agg_c1 / 4.0);   // beginning of first stage pi model
      stamp(Centries, 7, vic_c1 / 4.0);
      stamp(Gentries, 1, 2, 2/agg_r1);    // first stage pi resistance
      stamp(Gentries, 7, 8, 2/vic_r1);
      stamp(Centries, 2, agg_c1 / 4.0);   // end of first stage pi model
      stamp(Centries, 8, vic_c1 / 4.0);

      stamp(Centries, 2, agg_c1 / 4.0);   // beginning of second stage pi model
      stamp(Centries, 8, vic_c1 / 4.0);
      stamp(Gentries, 2, 3, 2/agg_r1);    // second stage pi resistance
      stamp(Gentries, 8, 9, 2/vic_r1);
      stamp(Centries, 3, agg_c1 / 4.0);   // end of second stage pi model
      stamp(Centries, 9, vic_c1 / 4.0);

      stamp(Centries, 3, 9, coup_c);      // coupling capacitance between traces

      stamp(Centries, 3, agg_c2 / 4.0);   // beginning of third stage pi model
      stamp(Centries, 9, vic_c2 / 4.0);
      stamp(Gentries, 3, 4, 2/agg_r2);    // third stage pi resistance
      stamp(Gentries, 9, 10, 2/vic_r2);
      stamp(Centries, 4, agg_c2 / 4.0);   // end of third stage pi model
      stamp(Centries, 10, vic_c2 / 4.0);

      stamp(Centries, 4, agg_c2 / 4.0);   // beginning of fourth stage pi model
      stamp(Centries, 10, vic_c2 / 4.0);
      stamp(Gentries, 4, 5, 2/agg_r2);    // fourth stage pi resistance
      stamp(Gentries, 10, 11, 2/vic_r2);
      stamp(Centries, 5, agg_c2 / 4.0);   // end of fourth stage pi model

      stamp(Centries, 5, agg_cl);
      // The victim receiver is the third port of the reduced network
      // we cannot put this through the reduction, because we model the ports as voltage
      // sources, which neutralizes capacitors.  Instead we will add this at the end.
      // stamp(Centries, 11, vic_cl);
      // stamp(Centries, 11, vic_c2 / 4.0);   // end of fourth stage pi model, victim side

      // add an additional pair of equations for the independent sources
      // aggressor and victim driver source currents will be nodes 12 and 13
      stamp_i(Gentries, 0, 12);
      stamp_i(Gentries, 6, 13);
      stamp_i(Gentries, 11, 14);           // victim receiver "driver" (we won't use it) current

      Eigen::SparseMatrix<double> C(15, 15); C.setFromTriplets(Centries.begin(), Centries.end());
      Eigen::SparseMatrix<double> G(15, 15); G.setFromTriplets(Gentries.begin(), Gentries.end());

      // Apply PRIMA reduction to the input model

      // Step 1: create B and L (input and output) matrices
      // In the PRIMA paper all ports are treated as inputs
      // Some postprocessing by the caller will fix this
      SparseMatrix<double> B(C.rows(), N);

      // insert a negative identity matrix in the last <N> rows
      B.reserve(N);     // 1 non-zero element per column
      size_t firstrow = C.rows() - N;
      for (size_t col = 0; col < N; ++col) {
        B.insert(firstrow + col, col) = -1;
      }

      SparseMatrix<double> L = -B;   // inputs are also outputs

      auto Xfinal = EDASkel::analysis::Prima(C, G, B, L, q);

      // Compute the C and G matrices in the new state variables using X
      Matrix<double, q, q> Cprime = Xfinal.transpose() * C * Xfinal;  // TODO: look at SVD here
      Matrix<double, q, q> Gprime = Xfinal.transpose() * G * Xfinal;
      Matrix<double, q, N> Bprime = Xfinal.transpose() * B;
      Matrix<double, q, N> Lprime = Xfinal.transpose() * L;

      // evaluate results

      // Compare block moments between original and reduced model
      // Prima claims to produce the same moments up to floor(q/N)
      SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > G_QR(G);
      assert(G_QR.info() == Success);
      // warning: in Pre-Eigen-3.2 debug builds this triggers an assertion failure due to an Eigen bug
      // see https://forum.kde.org/viewtopic.php?f=74&t=117474
      SparseMatrix<double> R = G_QR.solve(B);
      assert(G_QR.info() == Success);
      SparseMatrix<double> A = G_QR.solve(C);

      auto GprimeQR = Gprime.fullPivHouseholderQr();
      Matrix<double, q, q> Aprime = GprimeQR.solve(Cprime);
      Matrix<double, q, N> Rprime = GprimeQR.solve(Bprime);
      SparseMatrix<double> AtotheI(15, 15); AtotheI.setIdentity();

      Matrix<double, q, q> AprimetotheI = Matrix<double, q, q>::Identity();
      for (size_t i = 0; i < (q/N); ++i) {
        std::cerr << "moment " << i << " of original model is:" << std::endl;
        std::cerr << L.transpose() * AtotheI * R << std::endl;
        std::cerr << "moment " << i << " of reduced model is:" << std::endl;
        std::cerr << Lprime.transpose() * AprimetotheI * Rprime << std::endl;
        AtotheI = A * AtotheI;
        AprimetotheI = Aprime * AprimetotheI;
      }
      
      // Also compare eigenvalues, which are evidently the reciprocals of the poles
      // They are not guaranteed to be the same but should be "similar" (?)
      std::cerr << "eigenvalues of original model are:\n" << EigenSolver<MatrixXd>(A).eigenvalues() << std::endl;
      std::cerr << "eigenvalues of reduced model are:\n" << EigenSolver<decltype(Aprime)>(Aprime).eigenvalues() << std::endl;

      // PRIMA done.  Next step: Construct "Direct Stamp" reduced realization
      // The system equation (see "(50)" in the paper) hooks up the ports of
      // the reduced network to the remainder of the system for simulation
      // Fortunately there is no remainder for us as we have presented
      // the entire thing to PRIMA.  However, we will still need to turn the
      // victim receiver load (treated by PRIMA like all ports: as a voltage source)
      // into a plain circuit node so nothing flows out of the system and we
      // can simply observe the output.

      // state variables in our direct system are:
      // 0, 1, 2: node voltages for the two drivers and the victim receiver
      //          (two are directly driven while one will "float")
      // 3, 4: source currents for the two drivers
      // 5: "source current" for the last "driver" - always 0 as this is a receiver
      // 6.. : reduced system internal state variables (of quantity q)

      // equations for our direct system will be:
      // 2 equations to connect nodes 0 and 1 to the driver voltages
      // 1 equation connecting the victim receiver current to a small load
      // 3 equations connecting the reduced state to the current "outputs"
      // q equations for the normal reduced state evolution (will use state vars 0-2
      // for the input excitation)
      typedef Matrix<double, 6+q, 6+q> direct_eqn_matrix_t;
      direct_eqn_matrix_t Gdirect = direct_eqn_matrix_t::Zero();   // multiplies state vector
      direct_eqn_matrix_t Cdirect = direct_eqn_matrix_t::Zero();   // multiplies 1st derivative of state
      Matrix<double, 6+q, 2> Bdirect;
      Bdirect.block<2, 2>(0, 0) = -Matrix<double, 2, 2>::Identity();  // connects sources to nodes
      Matrix<double, 6+q, 1> Ldirect = Matrix<double, 6+q, 1>::Zero();
      Ldirect(2, 0) = 1;              // prima port 2 voltage is the sole output

      // Generate equations of the form Cdirect*(d(state)/dt) = -Gdirect*state + Bdirect*input;

      // hook up independent sources in first two rows: connects inputs to Vagg/Vvic via Bdirect
      Gdirect.block<2, 2>(0, 0) = -Matrix<double, 2, 2>::Identity(); 

      // Add back the victim receiver load capacitor and the final capacitor of the trace
      stamp(Cdirect, 2, vic_cl);
      stamp(Cdirect, 2, vic_c2 / 4.0);
      Gdirect(2, 5) = 1;              // connect port current to the caps

      // connect reduced state to currents (line 3 of "(50)")
      Gdirect.block<3, 3>(3, 3) = Matrix<double, 3, 3>::Identity();  // extract source currents
      Gdirect.block<3, q>(3, 6) = -Lprime.transpose();

      // connect reduced system  (line 4 of "(50)")

      auto Gprime_inv_x_Bprime = GprimeQR.solve(Bprime);
      auto Gprime_inv_x_Cprime = GprimeQR.solve(Cprime);

      Gdirect.block<q, 3>(6, 0) = Gprime_inv_x_Bprime;
      Gdirect.block<q, q>(6, 6) = Matrix<double, q, q>::Identity();
      Cdirect.block<q, q>(6, 6) = Gprime_inv_x_Cprime;
      
      // Direct Stamp realization constructed

      // reduce matrices to standard form so we can perform numeric integration
      MatrixXd Gred, Cred;
      Matrix<double, Dynamic, 2> Bred;
      Matrix<double, Dynamic, 1> Lred;
      std::tie(Gred, Cred, Bred, Lred) = regularize_su(Gdirect, Cdirect, Bdirect, Ldirect);

      // Solve to produce a pair of matrices we can use to calculate dX/dt
      auto CredQR = Cred.fullPivHouseholderQr();
      coeff_ = CredQR.solve(-1.0 * Gred);  // state evolution
      input_ = CredQR.solve(Bred);         // input -> state
      output_ = Lred.transpose();          // state -> output
  }

  void operator() (const state_type x, state_type& dxdt, double t) {
    // determine input voltages
    double Vagg = (agg_slew_ < 0) ? v_ : 0;   // initial value
    if (agg_slew_ == 0) {
      Vagg = 0.0;                             // quiescent
    } else {
      // find the "real" ramp starting point, since we use the center of the swing
      double transition_time = ( v_ / fabs(agg_slew_) ) / 2.0;
      double real_start = agg_start_ - transition_time;
      double real_end = agg_start_ + transition_time;

      if (t >= real_end) {
        Vagg = (agg_slew_ < 0) ? 0 : v_;      // straight to final value
      } else if (t <= real_start) {
        Vagg = (agg_slew_ > 0) ? 0 : v_;      // remain at initial value
      } else {
        Vagg += agg_slew_ * (t - real_start);  // proportional to time in ramp
      }
    }
        
    // same for the victim driver
    double Vvic = (vic_slew_ < 0) ? v_ : 0;   // initial value
    if (vic_slew_ == 0) {
      Vvic = 0.0;
    } else {
      // find the "real" ramp starting point, since we use the center of the swing
      double transition_time = ( v_ / fabs(vic_slew_) ) / 2.0;
      double real_start = vic_start_ - transition_time;
      double real_end = vic_start_ + transition_time;
      if (t >= real_end) {
        Vvic = (vic_slew_ < 0) ? 0 : v_;      // straight to final value
      } else if (t <= real_start) {
        Vvic = (vic_slew_ > 0) ? 0 : v_;      // remain at initial value
      } else {
        Vvic += vic_slew_ * (t - real_start);  // proportional to time in ramp
      }
    }
        
    // turn state vector into Eigen matrix
    Map<const Matrix<double, Dynamic, 1> > xvec(x.data(), statecount(), 1);

    // Input excitation:  the independent variables calculated above
    Matrix<double, 2, 1> u; u << Vagg, Vvic;

    Map<Matrix<double, Dynamic, 1> > result(dxdt.data(), statecount(), 1);
    result = coeff_ * xvec + input_ * u;              // sets dxdt via reference

  }

  Matrix<double, 1, q+1> output() const
  {
     // supply a matrix that can be used to extract the single output from the current state
     return output_;
  }

  // the actual number of state bits used in simulation depends on
  // the reduction results and on regularizing the circuit:
  size_t statecount() const {
    return coeff_.rows();
  }

};

// observer to record times and state values
struct push_back_state_and_time
{
    vector< state_type >& m_states;
    vector< double >& m_times;

    push_back_state_and_time( vector< state_type > &states , vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

// some measurement routines

double max_voltage(const vector<state_type>& history, size_t idx) {
  return (*max_element(history.begin(), history.end(),
		       [idx](state_type helt1, state_type helt2) {
			 return helt1[idx] < helt2[idx]; }))[idx];
}

enum SignalPairDirections { RiseRise, FallFall, RiseFall, FallRise };

// measure delay between two signals
double delay(const vector<double>& times,
	     const vector<state_type>& history,
	     double measurement_point,
	     size_t start_idx,
	     size_t stop_idx,
	     SignalPairDirections dirs) {
  // locate the first point at which the first (reference) signal passes the requested voltage
  // BOZO this would be a good place to do some error checking
  auto start_point = find_if(history.begin(), history.end(),
			     [start_idx, dirs, measurement_point](state_type helt) {
			       if ((dirs == RiseRise) || (dirs == RiseFall))
				 return helt[start_idx] >= measurement_point;
			       else
				 return helt[start_idx] <= measurement_point;
			     });
  double start_time = times[distance(history.begin(), start_point)];
  
  auto stop_point = find_if(history.begin(), history.end(),
			    [stop_idx, dirs, measurement_point](state_type helt) {
			      if ((dirs == RiseRise) || (dirs == RiseFall))
				 return helt[stop_idx] >= measurement_point;
			       else
				 return helt[stop_idx] <= measurement_point;
			     });
  double stop_time = times[distance(history.begin(), stop_point)];
  
  return stop_time - start_time;
}



int main() {
  // pick some seemingly sensible values
  double v = 1.0;
  double slew = v / 200e-12;   // 200ps rise/fall times
  double drvr_r = 100;
  double pi_r = 1000;          // resistance per segment
  double pi_c = 100e-15;       // 100fF
  double coupling_c = 100e-15;
  double rcvr_c = 20e-15;
  double drvr_start = 100e-12;

  signal_coupling ckt(pi_r, pi_c, pi_r, pi_c, rcvr_c, slew, drvr_r, drvr_start,
		      pi_r, pi_c, pi_r, pi_c, rcvr_c, 0.0, drvr_r, drvr_start,  // quiescent victim
		      coupling_c, v);

  // initial state: all low
  size_t statecount = ckt.statecount();
  state_type x(statecount, 0.0);

  vector<state_type> state_history;
  vector<state_type> outputs;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 1000e-12, 1e-12,
                     push_back_state_and_time( state_history, times ) );

  // generate observable (output) values by applying output transform matrix to each state
  auto const output_translator = ckt.output();
  transform(state_history.begin(), state_history.end(), back_inserter(outputs),
            [statecount, &output_translator]
            (const state_type& state) -> state_type {
              // access vector state data through an Eigen map so we can multiply
              Map<const Matrix<double, Dynamic, 1> > statevec(state.data(), statecount, 1);
              state_type ovec(1);   // result
              Map<Matrix<double, 1, 1> > outputvec(ovec.data(), 1, 1);   // also needs a Map
              outputvec = output_translator * statevec;  // multiply into vector via Map
              return ovec;
            });

  for (size_t i = 0; i < times.size(); ++i) {
    // format for gnuplot.  We are remembering the time and vvic only
     cout << times[i] << " " << outputs[i][0] << endl;
  }

  // find the highest voltage on the victim (which is supposed to be low)
  cerr << "max victim excursion is: " << max_voltage(outputs, 0) << endl;

  return 0;

}

const size_t signal_coupling::q;
