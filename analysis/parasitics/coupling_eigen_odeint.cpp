// Test of odeint using intersignal coupling
// Now with the use of Eigen to separate the equations
// Jeff Trull 2013-05-17

#include <vector>
using namespace std;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric;
#include <Eigen/Dense>
using namespace Eigen;

// Functions for implementing MNA with an Eigen matrix
template<typename M, typename Float>
void stamp(M& matrix, std::size_t i, std::size_t j, Float g)
{
   // General stamp: conductance at [i,i] and [j,j],
   // -conductance at [i,j] and [j,i], summed with existing values

   matrix(i, i) += g;
   matrix(j, j) += g;
   matrix(i, j) -= g;
   matrix(j, i) -= g;

}

// for when the other end of the device is at GND
template<typename M, typename Float>
void stamp(M& matrix, std::size_t i, Float g)
{
   matrix(i, i) += g;
}

// for voltage sources (inputs)
template<typename M>
void stamp_i(M& matrix, std::size_t vnodeno, std::size_t istateno)
{
   // just basically marks the connection between the inductor (or voltage source)
   // and the voltage, because unlike capacitance, both are state variables.
   // Doing this brings in the V = LdI/dt equations:

   matrix(vnodeno, istateno) = 1;   // current is taken *into* inductor or vsource
   matrix(istateno, vnodeno) = -1;
}

typedef vector<double> state_type;

struct signal_coupling {
  typedef Matrix<double, 6, 6> Matrix6d;
  Matrix6d coeff_;
  Matrix<double, 6, 2> input_;

  Matrix<double, 6, 3> Lred_;   // (regularized/reduced) state variable to output mapping
  Matrix<double, 3, 2> Dred_;   // regularized input to output mapping

  double agg_r1_, agg_c1_;   // aggressor first stage pi model (prior to coupling point)
  double agg_r2_, agg_c2_;   // aggressor second stage pi model (after coupling point)
  double agg_cl_;            // aggressor final load cap
  double agg_slew_;          // aggressor slew rate (V/s)
  double agg_imp_;           // aggressor driver impedance (placed after voltage source)
  double agg_start_;         // aggressor driver start time (the *center* of the ramp!)

  double vic_r1_, vic_c1_;   // victim first stage pi model (prior to coupling point)
  double vic_r2_, vic_c2_;   // victim second stage pi model (after coupling point)
  double vic_cl_;            // victim final load cap
  double vic_slew_;          // victim input slew rate (V/s)
  double vic_imp_;           // victim driver impedance
  double vic_start_;         // victim driver start time

  double coup_c_;            // coupling capacitance placed at central point
  double v_;                 // power supply (and thus max) voltage

  signal_coupling(double agg_r1, double agg_c1, double agg_r2, double agg_c2,
		  double agg_cl, double agg_slew, double agg_imp, double agg_start,
		  double vic_r1, double vic_c1, double vic_r2, double vic_c2,
		  double vic_cl, double vic_slew, double vic_imp, double vic_start,
		  double coup_c, double v) :
    agg_r1_(agg_r1), agg_c1_(agg_c1), agg_r2_(agg_r2), agg_c2_(agg_c2),
    agg_cl_(agg_cl), agg_slew_(agg_slew), agg_imp_(agg_imp), agg_start_(agg_start),
    vic_r1_(vic_r1), vic_c1_(vic_c1), vic_r2_(vic_r2), vic_c2_(vic_c2),
    vic_cl_(vic_cl), vic_slew_(vic_slew), vic_imp_(vic_imp), vic_start_(vic_start),
    coup_c_(coup_c), v_(v) {
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
      typedef Matrix<double, 10, 10> Matrix10d;
      Matrix10d C = Matrix10d::Zero(), G = Matrix10d::Zero();

      stamp(G, 0, 1, 1/agg_imp_);   // driver impedances
      stamp(G, 4, 5, 1/vic_imp_);

      stamp(C, 1, agg_c1_ / 2.0);   // beginning of first stage pi model
      stamp(C, 5, vic_c1_ / 2.0);

      stamp(G, 1, 2, 1/agg_r1_);    // first stage pi resistance
      stamp(G, 5, 6, 1/vic_r1_);
   
      stamp(C, 2, agg_c1_ / 2.0);   // end of first stage pi model
      stamp(C, 6, vic_c1_ / 2.0);

      stamp(C, 2, 6, coup_c_);      // coupling capacitance between traces

      stamp(C, 2, agg_c2_ / 2.0);   // beginning of second stage pi model
      stamp(C, 6, vic_c2_ / 2.0);

      stamp(G, 2, 3, 1/agg_r2_);    // second stage pi resistance
      stamp(G, 6, 7, 1/vic_r2_);
   
      stamp(C, 3, agg_c2_ / 2.0);   // end of second stage pi model
      stamp(C, 7, vic_c2_ / 2.0);

      stamp(C, 3, agg_cl_);
      stamp(C, 7, vic_cl_);

      // add an additional pair of equations for the independent sources
      // aggressor and victim driver source currents will be nodes 8 and 9
      stamp_i(G, 0, 8);
      stamp_i(G, 4, 9);

      // Model inputs and outputs with two matrices so that:
      // C*dX/dt = -G*X + B*u(t) and
      // Y = L.transpose() * X 

      // We have two inputs, so u(t) is 2x1 and B is 10x2
      Matrix<double, 10, 2> B; B << 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , -1, 0     // insert Vagg = V0
                                  , 0, -1 ;   // insert Vvic = V4

      // Similarly, two outputs (voltages at the receivers, not source currents!)
      Matrix<double, 10, 3> L; L << 0, 0, 0
                                  , 1, 0, 0    // extract V1 (aggressor driver output)
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 1, 0    // extract v6 (victim coupling node)
                                  , 0, 0, 1    // extract v7 (victim rcvr)
                                  , 0, 0, 0
                                  , 0, 0, 0 ;
                                  
      // 10x10 matrix from MNA cannot be directly integrated because some state
      // variable derivatives are not present.  This process substitutes other
      // state variables and in the process reduces the state size by 4 (by
      // eliminating the input voltage and current).
      // See comments in rlc_eigen_odeint.cpp

      // 1. Create a permutation of the state vector variables so those with
      //    non-zero rows in C are clustered at the top

      // Use Eigen reductions to find zero rows
      auto zero_rows = (C.array() == 0.0).rowwise().all();   // per row "all zeros"

      Matrix10d permut = Matrix10d::Identity();   // null permutation to start
      std::size_t i, j;
      for (i = 0, j=9; i < j;) {
        // loop invariant: rows > j are all zero; rows < i are not
        while ((i < 10) && !zero_rows(i)) ++i;
        while ((j > 0) && zero_rows(j)) --j;
        if (i < j) {
          // exchange rows i and j via the permutation vector
          permut(i, i) = 0; permut(j, j) = 0;
          permut(i, j) = 1; permut(j, i) = 1;
          ++i; --j;
        }
      }

      // 2. Apply permutation to MNA matrices
      Matrix10d Cprime = permut * C * permut;       // permute rows and columns
      Matrix10d Gprime = permut * G * permut;
      Matrix<double, 10, 2> Bprime = permut * B;    // permute only rows
      Matrix<double, 10, 3> Lprime = permut * L;
      
      // 3. Produce reduced equations following Su (Proc. 15th ASP-DAC, 2002)
      std::size_t zero_count = zero_rows.count();
      std::size_t nonzero_count = 10 - zero_count;

      auto G11 = Gprime.topLeftCorner(nonzero_count, nonzero_count);
      auto G12 = Gprime.topRightCorner(nonzero_count, zero_count);
      auto G21 = Gprime.bottomLeftCorner(zero_count, nonzero_count);
      auto G22 = Gprime.bottomRightCorner(zero_count, zero_count);

      auto L1 = Lprime.topRows(nonzero_count);
      auto L2 = Lprime.bottomRows(zero_count);

      auto B1 = Bprime.topRows(nonzero_count);
      auto B2 = Bprime.bottomRows(zero_count);

      auto Cred = Cprime.topLeftCorner(nonzero_count, nonzero_count);
      auto G22inv = G22.inverse();   // this is the most expensive thing we will do
      auto Gred = G11 - G12 * G22inv * G21;
      // our L, per PRIMA, is the transpose of the one described in Su
      Lred_ = (L1.transpose() - L2.transpose() * G22inv * G21).transpose();  // simplify?
      auto Bred = B1 - G12 * G22inv * B2;
      // we had no "D" (direct input to output) originally but it can happen
      Dred_ = L2.transpose() * G22inv * B2;

      // 4. Solve to produce a pair of matrices we can use to calculation dX/dt
      coeff_ = Cred.ldlt().solve(-1.0 * Gred);  // state evolution
      input_ = Cred.ldlt().solve(Bred);         // input -> state
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
        
    Map<const Matrix<double, 6, 1> > xvec(x.data());  // turn state vector into Eigen matrix
    Matrix<double, 2, 1> u; u << Vagg, Vvic;          // input excitation

    Map<Matrix<double, 6, 1> > result(dxdt.data());
    result = coeff_ * xvec + input_ * u;              // sets dxdt via reference

  }

  // utility function for displaying data.  Applies internal matrices to state
  std::vector<double> state2output(state_type const& x) {
    std::vector<double> result(2);
    Map<const Matrix<double, 6, 1> > xvec(x.data());

    Map<Matrix<double, 3, 1> > ovec(result.data());
    ovec = Lred_.transpose() * xvec;  // in theory the input could be involved via Dred_
    return result;
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
  state_type x(6, 0.0);

  vector<state_type> state_history;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 1000e-12, 1e-12,
                     push_back_state_and_time( state_history, times ) );

  for (size_t i = 0; i < times.size(); ++i) {
    // format for gnuplot.  look at driver output, victim coupling node, and victim receiver node
    auto output = ckt.state2output(state_history[i]);
    cout << times[i] << " " << output[0] << " " << output[1] << " " << output[2] << endl;
  }

  // find the highest voltage on the victim (which is supposed to be low)
  cerr << "max victim excursion is: " << max_voltage(state_history, 0) << endl;

  cerr << "driver delay is: " << delay(times, state_history, v/2.0, 1, 3, RiseRise) << endl;

  return 0;

}
