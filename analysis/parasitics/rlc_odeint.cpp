// Test of odeint using an RLC tank circuit
// Jeff Trull 2012-05-14

#include <vector>
using namespace std;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric;

using state_type = vector<double>;

struct rlc_tank {
  double r_, l_, c_;
  double vin_;

  rlc_tank(double r, double l, double c, double vin) : r_(r), l_(l), c_(c), vin_(vin) {}
  void operator() (const state_type x, state_type& dxdt, double t) {
    // based on manual algebra the right equation for series R, parallel LC
    // state variables are:
    // 0 - current through inductor
    // 1 - output voltage
    // 2 - input voltage
    dxdt[0] = x[1] / l_;                         // from V = L*dI/dt
    dxdt[1] = ( (x[2] - x[1])/r_ - x[0] ) / c_;  // from KCL and I = C*dV/dt
    // for input voltage, model an impulse at time zero
    if ((t < 0.0) || (x[2] >= vin_)) {
      dxdt[2] = 0.0;
    } else {
      dxdt[2] = numeric_limits<double>::max();
    }

    // this didn't work (differentiation eliminated necessary voltage source)
    // this strategy also tried to handle the step at time zero by supplying
    // 0 for the initial state of the output.  That may still work and avoid the
    // "impulse" hack
    /*
    // at output is:
    // V'' = -1/RC * (R/L * V + V')
    // using the state assignment x[0] = V, x[1] = V' we get:
    dxdt[0] = x[1];          // by definition
    // KCL on the output node produces this one
    //    dxdt[1] = -1.0/(r_*c_) * ( (r_/l_)*x[0] + x[1] );
    // introduce another variable for the excitation voltage
    */

  }
};

// observer to record times and state values
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};


int main() {
  // pick some seemingly sensible values
  double v = 3.3;
  double r = 100;
  double l = 20.0e-6;   // 20 uH
  double c = 20.0e-9;   // 20 nF
  rlc_tank ckt(r, l, c, v);


  // initial state: all voltages and currents 0
  state_type x(3, 0.0);

  vector<state_type> state_history;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 10e-6, 1e-6,
                     push_back_state_and_time( state_history, times ) );

  for (size_t i = 0; i < times.size(); ++i) {
    // format for gnuplot.  Only vin and vout are interesting
    cout << times[i] << " " << state_history[i][1] << " " << state_history[i][2] << endl;
  }

  return 0;

}
