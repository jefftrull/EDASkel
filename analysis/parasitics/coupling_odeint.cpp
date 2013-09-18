// Test of odeint using intersignal coupling
// Jeff Trull 2012-05-16

#include <vector>
using namespace std;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric;

typedef vector<double> state_type;

struct signal_coupling {
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
    coup_c_(coup_c), v_(v) {}
				   
  void operator() (const state_type x, state_type& dxdt, double t) {
    // state mapping (all voltages):

    // 0 - aggressor voltage source
    dxdt[0] = 0.0;
    // slew of 0.0 means constant
    if (agg_slew_ != 0.0) {
      // find the "real" ramp starting point, since we use the center of the swing
      double real_start = agg_start_ - ( ( v_ / fabs(agg_slew_) ) / 2.0 );
      if ((t >= real_start) &&
	  (((agg_slew_ > 0.0) && (x[0] < v_)) ||
	   ((agg_slew_ < 0.0) && (x[0] > 0.0)))) {
	  // we are past the starting point but haven't reached our final voltage
	  dxdt[0] = agg_slew_;
	}
    }

    // 1 - aggressor driver output (past the impedance model resistor)
    double cdrv_agg = agg_c1_ / 2.0;
    dxdt[1] = ( ( (x[0] - x[1]) / agg_imp_ ) -              // KCL for capacitor current
		( (x[1] - x[2]) / agg_r1_ ) ) / cdrv_agg;  // divided by capacitance

    // 2 - aggressor coupling node - tricky because the currents on the coupled nodes
    //     are interrelated.  Fortunately we have two equations (KCL) in two unknowns
    //     (dV/dt of each side) and we can solve for each dV/dt

    //     net resistor current into node (all other currents depend on dvdt)
    double ires_agg = ( (x[1] - x[2]) / agg_r1_ ) - ( (x[2] - x[3]) / agg_r2_);
    //     we need the victim one too
    double ires_vic = ( (x[5] - x[6]) / vic_r1_ ) - ( (x[6] - x[7]) / vic_r2_);
    //     ground-lumped (non-coupling) capacitive load
    double clmp_agg = ( agg_c1_ + agg_c2_ ) / 2.0;
    double clmp_vic = ( vic_c1_ + vic_c2_ ) / 2.0;
    //     coupling ratios (relates lumped and coupling)
    double coupr_agg = (clmp_agg + coup_c_) / clmp_agg;
    double coupr_vic = (clmp_vic + coup_c_) / clmp_vic;
    //     now the equation for dVa/dt - manually derived
    dxdt[2] = ( ( coupr_vic * (ires_vic + ires_agg) - ires_vic ) /
		( coupr_vic * clmp_agg + coup_c_ ) );
    
    // 3 - aggressor receiver node
    dxdt[3] = ( ( ( x[2] - x[3] ) / agg_r2_ ) /       // current into receiver
		( (agg_c2_ / 2.0) + agg_cl_) );       // load at receiver

    // 4 - victim driving node
    dxdt[4] = 0.0;
    // slew of 0.0 means constant
    if (vic_slew_ != 0.0) {
      // find the "real" ramp starting point, since we use the center of the swing
      double real_start = vic_start_ - ( ( v_ / fabs(vic_slew_) ) / 2.0 );
      if ((t >= real_start) &&
	  (((vic_slew_ > 0.0) && (x[0] < v_)) ||
	   ((vic_slew_ < 0.0) && (x[0] > 0.0)))) {
	  // we are past the starting point but haven't reached our final voltage
	  dxdt[4] = vic_slew_;
	}
    }

    double real_vic_start = vic_start_ - ( ( v_ / fabs(vic_slew_) ) / 2.0 );
    dxdt[4] = 0.0;
    if (t >= real_vic_start) {
      if (((vic_slew_ > 0.0) && (x[4] < v_)) ||
	  ((vic_slew_ < 0.0) && (x[4] > 0.0))) {
	dxdt[4] = vic_slew_;
      }
    }

    // 5 - victim driver output
    double cdrv_vic = vic_c1_ / 2.0;
    dxdt[5] = ( ( (x[4] - x[5]) / vic_imp_ ) -             // KCL for capacitor current
		( (x[5] - x[6]) / vic_r1_ ) ) / cdrv_vic;  // divided by capacitance

    // 6 - victim coupling node
    dxdt[6] = ( ( coupr_agg * (ires_agg + ires_vic) - ires_agg ) /
		( coupr_agg * clmp_vic + coup_c_ ) );

    // 7 - victim receiver node
    dxdt[7] = ( ( ( x[6] - x[7] ) / vic_r2_ ) /       // current into receiver
		( (vic_c2_ / 2.0) + vic_cl_) );       // load at receiver

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
  state_type x(8, 0.0);

  vector<state_type> state_history;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 1000e-12, 1e-12,
                     push_back_state_and_time( state_history, times ) );

  for (size_t i = 0; i < times.size(); ++i) {
    // format for gnuplot.  look at driver output, victim coupling node, and victim receiver node
    cout << times[i] << " " << state_history[i][1] << " " << state_history[i][6] << " " << state_history[i][7] << endl;
    //    cout << times[i] << " " << state_history[i][1] << " " << state_history[i][3] << endl;
  }

  // find the highest voltage on the victim (which is supposed to be low)
  cerr << "max victim excursion is: " << max_voltage(state_history, 7) << endl;

  cerr << "driver delay is: " << delay(times, state_history, v/2.0, 1, 3, RiseRise) << endl;

  return 0;

}
