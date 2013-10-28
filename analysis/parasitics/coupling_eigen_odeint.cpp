// Test of odeint using intersignal coupling
// Now with the use of Eigen to separate the equations
// Jeff Trull 2013-05-17

#include <vector>
using namespace std;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric;
#include <Eigen/Dense>
using namespace Eigen;

#include "analysis/mna.hpp"

typedef vector<double> state_type;

struct signal_coupling {
  static const size_t states = 10;

  MatrixXd coeff_;
  MatrixXd input_;

  MatrixXd Lred_;            // (regularized/reduced) state variable to output mapping
  MatrixXd Dred_;            // regularized input to output mapping

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
      typedef Matrix<double, 14, 14> Matrix14d;
      Matrix14d C = Matrix14d::Zero(), G = Matrix14d::Zero();

      using namespace EDASkel::analysis::mna;
      stamp(G, 0, 1, 1/agg_imp_);   // driver impedances
      stamp(G, 6, 7, 1/vic_imp_);

      stamp(C, 1, agg_c1_ / 4.0);   // beginning of first stage pi model
      stamp(C, 7, vic_c1_ / 4.0);
      stamp(G, 1, 2, 2/agg_r1_);    // first stage pi resistance
      stamp(G, 7, 8, 2/vic_r1_);
      stamp(C, 2, agg_c1_ / 4.0);   // end of first stage pi model
      stamp(C, 8, vic_c1_ / 4.0);

      stamp(C, 2, agg_c1_ / 4.0);   // beginning of second stage pi model
      stamp(C, 8, vic_c1_ / 4.0);
      stamp(G, 2, 3, 2/agg_r1_);    // second stage pi resistance
      stamp(G, 8, 9, 2/vic_r1_);
      stamp(C, 3, agg_c1_ / 4.0);   // end of second stage pi model
      stamp(C, 9, vic_c1_ / 4.0);

      stamp(C, 3, 9, coup_c_);      // coupling capacitance between traces

      stamp(C, 3, agg_c2_ / 4.0);   // beginning of third stage pi model
      stamp(C, 9, vic_c2_ / 4.0);
      stamp(G, 3, 4, 2/agg_r2_);    // third stage pi resistance
      stamp(G, 9, 10, 2/vic_r2_);
      stamp(C, 4, agg_c2_ / 4.0);   // end of third stage pi model
      stamp(C, 10, vic_c2_ / 4.0);

      stamp(C, 4, agg_c2_ / 4.0);   // beginning of fourth stage pi model
      stamp(C, 10, vic_c2_ / 4.0);
      stamp(G, 4, 5, 2/agg_r2_);    // fourth stage pi resistance
      stamp(G, 10, 11, 2/vic_r2_);
      stamp(C, 5, agg_c2_ / 4.0);   // end of fourth stage pi model
      stamp(C, 11, vic_c2_ / 4.0);

      stamp(C, 5, agg_cl_);
      stamp(C, 11, vic_cl_);

      // add an additional pair of equations for the independent sources
      // aggressor and victim driver source currents will be nodes 12 and 13
      stamp_i(G, 0, 12);
      stamp_i(G, 6, 13);

      // Model inputs and outputs with two matrices so that:
      // C*dX/dt = -G*X + B*u(t) and
      // Y = L.transpose() * X 

      // We have two inputs, so u(t) is 2x1 and B is 14x2
      Matrix<double, 14, 2> B; B << 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , -1, 0     // insert Vagg = V0
                                  , 0, -1 ;   // insert Vvic = V6

      // Three outputs, for viewing and performing measurements
      Matrix<double, 14, 3> L; L << 0, 0, 0
                                  , 1, 0, 0    // extract V1 (aggressor driver output)
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 1, 0    // extract v7 (victim coupling node)
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 1    // extract v11 (victim rcvr)
                                  , 0, 0, 0
                                  , 0, 0, 0 ;
                                  
      // 14x14 matrix from MNA cannot be directly integrated because some state
      // variable derivatives are not present.  This process substitutes other
      // state variables and in the process reduces the state size by 4 (by
      // eliminating the input voltage and current).
      // See comments in rlc_eigen_odeint.cpp

      MatrixXd Cred, Gred;
      Matrix<double, Dynamic, 2> Bred;
      std::tie(Gred, Cred, Bred, Lred_) = regularize(G, C, B, L);

      // Solve to produce a pair of matrices we can use to calculate dX/dt
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
        
    Map<const Matrix<double, states, 1> > xvec(x.data());  // turn state vector into Eigen matrix
    Matrix<double, 2, 1> u; u << Vagg, Vvic;          // input excitation

    Map<Matrix<double, states, 1> > result(dxdt.data());
    result = coeff_ * xvec + input_ * u;              // sets dxdt via reference

  }

  // utility function for displaying data.  Applies internal matrices to state
  std::vector<double> state2output(state_type const& x) {
    std::vector<double> result(2);
    Map<const Matrix<double, states, 1> > xvec(x.data());

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
  state_type x(signal_coupling::states, 0.0);

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

  cerr << "driver delay is: " << delay(times, state_history, v/2.0, 1, 5, RiseRise) << endl;

  return 0;

}
