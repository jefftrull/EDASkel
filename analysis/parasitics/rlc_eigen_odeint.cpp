// Test of odeint+Eigen(MNA) using an RLC tank circuit
// Jeff Trull 2013-05-17

#include <vector>
using namespace std;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric;
#include <Eigen/Dense>
using namespace Eigen;

// Functions for implementing MNA with an Eigen matrix
// general stamp for conductance
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

// conductance stamp for when the other end of the device is at GND
template<typename M, typename Float>
void stamp(M& matrix, std::size_t i, Float g)
{
   matrix(i, i) += g;
}

// for inductor currents
template<typename M>
void stamp_i(M& matrix, std::size_t vnodeno, std::size_t istateno)
{
   // just basically marks the connection between the inductor (or voltage source)
   // and the voltage, because unlike capacitance, both are state variables.
   // Doing this brings in the V = LdI/dt equations:

   matrix(vnodeno, istateno) = 1;   // current is taken *into* inductor or vsource
   matrix(istateno, vnodeno) = 1;
}

typedef vector<double> state_type;

struct rlc_tank {
  double r_, l_, c_;
  double vin_;
  Matrix4f coeff_;

  rlc_tank(double r, double l, double c, double vin) : r_(r), l_(l), c_(c), vin_(vin) {
    // state variables are:
    // 0 - input voltage
    // 1 - output voltage
    // 2 - current from input voltage source
    // 3 - current through inductor

    // Create matrices implementing the equation G*X + C*dX/dt = u(t)
    // where u(t) is the independent sources

    Matrix4f C = Matrix4f::Zero(), G = Matrix4f::Zero();

    stamp(G, 0, 1, 1.0/r_);
    stamp(C, 1, c_);
    stamp_i(G, 0, 2);   // connect input current to voltage source
    // for the inductor, two operations:
    stamp(C, 3, -l_);   // set derivative coefficient (minus, because same side as voltage)
    stamp_i(G, 1, 3);   // connect inductor current derivative to associated voltage

    // Rewrite equation as C*dX/dt = -G*X + u(t)
    // Solve to produce a single matrix with the equations for each node
    coeff_ = C.fullPivLu().solve(-1.0 * G);
     
  }

  void operator() (const state_type x, state_type& dxdt, double t) {

    // for input voltage, model a step function at time zero
    if ((t < 0.0) || (x[0] >= vin_)) {
      dxdt[0] = 0.0;
    } else {
      dxdt[0] = numeric_limits<double>::max();
    }

    // All other node voltages are determined by odeint through our equations:
    for (std::size_t nodeno = 1; nodeno <= 3; ++nodeno)
    {
      // BOZO std::accumulate with lambda or something
      dxdt[nodeno] = 0.f;
      for (std::size_t stateno = 0; stateno <= 3; ++stateno)
      {
         dxdt[nodeno] += coeff_(nodeno, stateno) * x[stateno];
      }
    }
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
  state_type x({0.0, 0.0, 0.0, 0.0});

  vector<state_type> state_history;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 10e-6, 1e-6,
                     push_back_state_and_time( state_history, times ) );

  for (size_t i = 0; i < times.size(); ++i) {
    // format for gnuplot.  Only vin and vout are interesting
    cout << times[i] << " " << state_history[i][1] << " " << state_history[i][0] << endl;
  }

  return 0;

}
