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
   matrix(istateno, vnodeno) = -1;
}

typedef vector<double> state_type;

struct rlc_tank {
  double r_, l_, c_;
  double vin_;
  MatrixXd coeff_;
  MatrixXd input_;

  MatrixXd Lred_;   // (regularized/reduced) state variable to output mapping
  MatrixXd Dred_;   // regularized input to output mapping

  rlc_tank(double r, double l, double c, double vin) : r_(r), l_(l), c_(c), vin_(vin) {
    // state variables are:
    // 0 - input voltage
    // 1 - output voltage
    // 2 - current from input voltage source
    // 3 - current through inductor

    // Create matrices implementing the equation G*X + C*dX/dt = u(t)
    // where u(t) is the independent sources

    Matrix4d C = Matrix4d::Zero(), G = Matrix4d::Zero();

    stamp(G, 0, 1, 1.0/r_);
    stamp(C, 1, c_);
    stamp_i(G, 0, 2);   // connect input current to voltage source
    // for the inductor, two operations:
    stamp(C, 3, l_);    // set derivative coefficient
    stamp_i(G, 1, 3);   // connect inductor current derivative to associated voltage

    // input application matrix
    Matrix<double, 4, 1> B; B << 0, 0, -1, 0;   // match Va up with V0

    // output observation matrix - we want to see input, output, and input current
    Matrix<double, 4, 3> L; L << 1, 0, 0
                               , 0, 1, 0
                               , 0, 0, -1       // extract input current as *out* of the source
                               , 0, 0, 0 ;

    // Rewrite equation as C*dX/dt = -G*X + B*u(t)
    // Further the variables we want to view are Y = L.transpose() * X

    // Now we have a set of differential equations, not all of which are in the required form
    // for ODEINT.  Specifically, some of the derivative terms have a coefficient of zero,
    // so they cannot be meaningfully integrated.  This problem is discussed in detail in:
    // Chen, "A Practical Regularization Technique for Modiﬁed Nodal Analysis...",
    // IEEE TCAD, July 2012 and my chosen solution is the simple one used in Su,
    // "Efﬁcient Approximate Balanced Truncation of General Large-Scale RLC Systems via Krylov Methods"
    // Proc. 15th ASP-DAC, 2002

    // Permute C to move the non-zero rows of C to the top while creating a "permutation" matrix
    // we can use to adjust G and B for the state variable reordering.

    // Use Eigen reductions to find zero rows
    auto zero_rows = (C.array() == 0.0).rowwise().all();   // per row "all zeros"

    Matrix4d permut = Matrix4d::Identity();   // null permutation to start
    std::size_t i, j;
    for (i = 0, j=3; i < j;) {
      // loop invariant: rows > j are all zero; rows < i are not
      while ((i < 4) && !zero_rows(i)) ++i;
      while ((j > 0) && zero_rows(j)) --j;
      if (i < j) {
        // exchange rows i and j via the permutation vector
        permut(i, i) = 0; permut(j, j) = 0;
        permut(i, j) = 1; permut(j, i) = 1;
        ++i; --j;
      }
    }

    // 2. Apply permutation to MNA matrices
    Matrix4d Cprime = permut * C * permut;       // permute rows and columns
    Matrix4d Gprime = permut * G * permut;
    Vector4d Bprime = permut * B;    // permute only rows
    Matrix<double, 4, 3> Lprime = permut * L;

    // now the first nonzero_count rows of Cprime, Gprime, and Bprime contain equations acceptable to
    // ODEINT, but the remaining rows do not.  We partition the state space into integrable state (X1)
    // and non-integrable (X2) and do the same with the matrices
    std::size_t zero_count = zero_rows.count();
    std::size_t nonzero_count = 4 - zero_count;

    auto G11 = Gprime.topLeftCorner(nonzero_count, nonzero_count);
    auto G12 = Gprime.topRightCorner(zero_count, zero_count);
    auto G21 = Gprime.bottomLeftCorner(zero_count, zero_count);
    auto G22 = Gprime.bottomRightCorner(zero_count, zero_count);

    auto L1 = Lprime.topRows(nonzero_count);
    auto L2 = Lprime.bottomRows(zero_count);

    auto B1 = Bprime.topRows(nonzero_count);
    auto B2 = Bprime.bottomRows(zero_count);

    // produce reduced equations following Su:
    auto Cred = Cprime.topLeftCorner(nonzero_count, nonzero_count);
    auto G22inv = G22.inverse();   // this is the most expensive thing we will do
    auto Gred = G11 - G12 * G22inv * G21;
    // our L, per PRIMA, is the transpose of the one described in Su
    Lred_ = (L1.transpose() - L2.transpose() * G22inv * G21).transpose();  // simplify?
    auto Bred = B1 - G12 * G22inv * B2;
    // we had no "D" (direct input to output) originally but it can happen, especially
    // since we directly control one state variable that we are mapping to an output
    Dred_ = L2.transpose() * G22inv * B2;

    // Solve to produce a single matrix with the equations for each node
    coeff_ = Cred.ldlt().solve(-1.0 * Gred);

    input_ = Cred.ldlt().solve(Bred);

  }

  void operator() (const state_type x, state_type& dxdt, double t) {

    // for input voltage, model a step function at time zero
    double Va = 0;
    if (t > 0.0) {
      Va = vin_;
    }
    Matrix<double, 1, 1> u; u << Va;

    // All other node voltages are determined by odeint through our equations:
    Map<const Matrix<double, 2, 1> > xvec(x.data());
    Map<Matrix<double, 2, 1> > result(dxdt.data());

    result = coeff_ * xvec + input_ * u;
  }

  // utility function for displaying data.  Applies internal matrices to state
  std::vector<double> state2output(state_type const& x,
                                   std::vector<double> const& in) {
    std::vector<double> result(3);
    Map<const Vector2d> xvec(x.data());
    Map<const Matrix<double, 1, 1> > invec(in.data());

    Map<Matrix<double, 3, 1> > ovec(result.data());
    ovec = Lred_.transpose() * xvec + Dred_ * invec;
    return result;
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
  state_type x(2, 0.0);

  vector<state_type> state_history;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 10e-6, 1e-6,
                     push_back_state_and_time( state_history, times ) );

  for (size_t i = 0; i < times.size(); ++i) {
    double Va = 0.0;
    if (times[i] > 0)
    {
      Va = v;
    }
    auto output = ckt.state2output(state_history[i], vector<double>(1, Va));

    // format for gnuplot.  Only vin and vout are interesting
    cout << times[i] << " " << output[1] << " " << output[0] << endl;
  }

  return 0;

}
