// Signal coupling case using SPRIM reduction
// Jeff Trull 2014-04-01

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
  static const size_t q = 6; // desired number of state variables in the reduced system
                             // 12 is the "natural" count (1 per circuit node)
  static const size_t N = 3; // port count (one per input/output: two drivers and the victim rcvr)

  // final system for simulation.  State variable count determined by regularization:
  Matrix<double, Dynamic, Dynamic> coeff_;  // reduced system state evolution
  Matrix<double, Dynamic, 2> input_;        // inputs (agg/vic) to reduced system state
  Matrix<double, 1, Dynamic> output_;       // reduced system state to chosen outputs
  Matrix<double, 1, 2> feedthrough_;        // direct contribution to output from inupt

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

      // Apply SPRIM reduction to the input model
      // This appears to be a simple matter of running PRIMA, then splitting the result into 2 or 3
      // (if there are inductors) matrices applied separately to different parts of the input

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

      // now do the postprocessing steps described in the SPRIM paper
      // We don't have any inductors, so only the first and last portions
      // (corresponding to components and the ports) are kept, i.e., states
      // 0-11 and 12-14.
      auto V1 = Xfinal.topRows(15 - N);
      auto V3 = Xfinal.bottomRows(N);
      // V3 generally lacks full column rank, so reduce its columns until this is true.
      auto V3QR = V3.fullPivHouseholderQr();
      auto V3QR_Q = V3QR.matrixQ();
      auto V3_fullrank = V3QR_Q.leftCols(V3QR.rank());
      std::size_t v3cols = V3_fullrank.cols();

      // G is split into D11 (12x12 corresponding to resistors) and Av
      // (Voltage source inputs), which is 12x3 to line up with D11 although
      // only the top 3x3 has data

      auto G11_prime = V1.transpose() * G.block(0, 0, 12, 12).selfadjointView<Upper>() * V1;  // D11
      auto G12_prime = V1.transpose() * G.block(0, 12, 12, N) * V3_fullrank;  // Va
      auto C11_prime = V1.transpose() * C.block(0, 0, 12, 12).selfadjointView<Upper>() * V1;  // E11

      // construct C/G in new state variables using SPRIM directions
      // TODO try using sparse matrices
      // The resulting matrices are square, with q rows/columns from the original
      // reduced system and another v3cols from the shifted and rank-reduced port block
      Matrix<double, Dynamic, Dynamic> Gprime, Cprime;
      std::size_t reducedstates = q + v3cols;
      Gprime = Matrix<double, Dynamic, Dynamic>::Zero(reducedstates, reducedstates);
      Gprime.topLeftCorner(q, q)         =  G11_prime;
      Gprime.topRightCorner(q, v3cols)   =  G12_prime;
      Gprime.bottomLeftCorner(v3cols, q) = -G12_prime.transpose();
      Cprime = Matrix<double, Dynamic, Dynamic>::Zero(reducedstates, reducedstates);
      Cprime.topLeftCorner(q, q)         =  C11_prime;
      
      Matrix<double, Dynamic, N> Bprime, Lprime;
      Bprime = Matrix<double, Dynamic, N>::Zero(reducedstates, N);
      Bprime.bottomRows(v3cols) = -V3_fullrank.transpose();
      Lprime               = -Bprime;

      // All of the above seems to actually work and produce the same result as multiplying

      // evaluate results

      // Compare block moments between original and reduced model
      // Prima claims to produce the same moments up to floor(q/N)
      // SPRIM claims to produce twice as many:
      Matrix<double, N, N> Ezero = Matrix<double, N, N>::Zero();  // no feedthrough
      auto moments_orig    = moments(Matrix<double, 15, 15>(G), Matrix<double, 15, 15>(C),
                                     Matrix<double, 15, N>(B), Matrix<double, 15, N>(L),
                                     Ezero, 2 * q/N);
      auto moments_reduced = moments(Gprime, Cprime, Bprime, Lprime,
                                     Ezero, 2 * q/N);

      for (size_t i = 0; i < (2 * q/N); ++i) {
        std::cerr << "moment " << i << " of original model is:" << std::endl;
        std::cerr << moments_orig[i] << std::endl;
        std::cerr << "moment " << i << " of reduced model is:" << std::endl;
        std::cerr << moments_reduced[i] << std::endl;
      }
      
      // Also compare eigenvalues, which are evidently the reciprocals of the poles
      // They are not guaranteed to be the same but should be "similar" (?)
      SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > G_QR(G);
      assert(G_QR.info() == Success);

      // This next step is what makes it pointless to use selfadj for C
      // according to ggael (and SparseSolverBase.h) a sparse RHS is converted into dense
      // temporary "panels" for solution
      SparseMatrix<double> A = G_QR.solve(C);
      assert(G_QR.info() == Success);
      MatrixXd Adense = A;
      std::cerr << "eigenvalues of original model are:\n" << EigenSolver<MatrixXd>(Adense).eigenvalues() << std::endl;
      auto GprimeQR = Gprime.fullPivHouseholderQr();
      Matrix<double, Dynamic, Dynamic> Aprime = GprimeQR.solve(Cprime);
      std::cerr << "eigenvalues of reduced model are:\n" << EigenSolver<decltype(Aprime)>(Aprime).eigenvalues() << std::endl;

      // Reduction complete.  Next step: Construct "Direct Stamp" reduced realization

      // state variables in our direct system are:
      // 0, 1, 2: node voltages for the two drivers and the victim receiver
      //          (two are directly driven while one will "float")
      // 3, 4: source currents for the two drivers
      // 5: "source current" for the last "driver" - always 0 as this is a receiver
      // 6.. : reduced system internal state variables (of quantity q+3)

      // equations for our direct system will be:
      // 2 equations to connect nodes 0 and 1 to the driver voltages
      // 1 equation connecting the victim receiver current to a small load
      // 3 equations connecting the reduced state to the current "outputs"
      // 2q equations for the normal reduced state evolution (will use state vars 0-2
      // for the input excitation)
      typedef Matrix<double, Dynamic, Dynamic> direct_eqn_matrix_t;
      direct_eqn_matrix_t Gdirect = direct_eqn_matrix_t::Zero(6+reducedstates, 6+reducedstates);   // multiplies state vector
      direct_eqn_matrix_t Cdirect = direct_eqn_matrix_t::Zero(6+reducedstates, 6+reducedstates);   // multiplies 1st derivative of state
      Matrix<double, Dynamic, 2> Bdirect = Matrix<double, Dynamic, 2>::Zero(6+reducedstates, 2);
      Bdirect.block<2, 2>(0, 0) = -Matrix<double, 2, 2>::Identity();  // connects sources to nodes
      Matrix<double, Dynamic, 1> Ldirect = Matrix<double, Dynamic, 1>::Zero(6+reducedstates, 1);
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
      Gdirect.block(3, 6, Lprime.cols(), Lprime.rows()) = -Lprime.transpose();

      // connect reduced system  (line 4 of "(50)")

      auto Gprime_inv_x_Bprime = GprimeQR.solve(Bprime);
      auto Gprime_inv_x_Cprime = GprimeQR.solve(Cprime);

      Gdirect.block(6, 0, Gprime_inv_x_Bprime.rows(), Gprime_inv_x_Bprime.cols()) = Gprime_inv_x_Bprime;
      Gdirect.block(6, 6, reducedstates, reducedstates) =
        Matrix<double, Dynamic, Dynamic>::Identity(reducedstates, reducedstates);
      Cdirect.block(6, 6, Gprime_inv_x_Cprime.rows(), Gprime_inv_x_Cprime.cols()) = Gprime_inv_x_Cprime;

      // calculate moments of direct system
      Matrix<double, 1, 2> Ezero12 = Matrix<double, 1, 2>::Zero();  // no feedthrough
      auto moments_direct = moments(Gdirect, Cdirect, Bdirect, Ldirect,
                                    Ezero12, q/N);
      for (size_t i = 0; i < (q/N); ++i) {
        std::cerr << "moment " << i << " of direct model is:" << std::endl;
        std::cerr << moments_direct[i] << std::endl;
      }

      // Direct Stamp realization constructed

      // reduce matrices to standard form so we can perform numeric integration
      MatrixXd Gred, Cred;
      Matrix<double, Dynamic, 2> Bred;
      Matrix<double, Dynamic, 1> Lred;
      Matrix<double, 1, 2> E;
      // We cannot use Su regularization because it requires the number of zero rows and
      // columns to be the same
      std::tie(Gred, Cred, Bred, Lred, E) = regularize_natarajan<2, 1, Dynamic>(Gdirect, Cdirect, Bdirect, Ldirect);

      // calculate moments of regularized direct system
      auto moments_reg = moments(Gred, Cred, Bred, Lred, E, q/N);
      for (size_t i = 0; i < (q/N); ++i) {
        std::cerr << "moment " << i << " of regularized model is:" << std::endl;
        std::cerr << moments_reg[i] << std::endl;
      }

      // Solve to produce a pair of matrices we can use to calculate dX/dt

      assert(!isSingular(Cred));
      auto CredQR = Cred.fullPivHouseholderQr();

      coeff_ = CredQR.solve(-1.0 * Gred);  // state evolution
      input_ = CredQR.solve(Bred);         // input -> state
      output_ = Lred.transpose();          // state -> output
      feedthrough_ = E;                    // input -> output

  }

  Matrix<double, 2, 1> inputs(double t) {
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
        
    // Input excitation:  the independent variables calculated above
    Matrix<double, 2, 1> u; u << Vagg, Vvic;

    return u;
  }

  void operator() (const state_type x, state_type& dxdt, double t) {
    // turn state vector into Eigen matrix
    Map<const Matrix<double, Dynamic, 1> > xvec(x.data(), statecount(), 1);

    Map<Matrix<double, Dynamic, 1> > result(dxdt.data(), statecount(), 1);
    result = coeff_ * xvec + input_ * inputs(t);           // sets dxdt via reference

  }

  Matrix<double, 1, Dynamic> output() const
  {
     // supply a matrix that can be used to extract the single output from the current state
     return output_;
  }

  Matrix<double, 1, 2> feedthrough() const
  {
      return feedthrough_;
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
  auto const state_to_output = ckt.output();
  auto const input_to_output = ckt.feedthrough();
  transform(state_history.begin(), state_history.end(),
            times.begin(), back_inserter(outputs),
            [statecount, &state_to_output, &input_to_output, &ckt]
            (const state_type& state, double t) -> state_type {
              // access vector state data through an Eigen map so we can multiply
              Map<const Matrix<double, Dynamic, 1> > statevec(state.data(), statecount, 1);
              state_type ovec(1);   // result
              Map<Matrix<double, 1, 1> > outputvec(ovec.data(), 1, 1);   // also needs a Map
              outputvec = state_to_output * statevec  // multiply into vector via Map
                        + input_to_output * ckt.inputs(t);
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
