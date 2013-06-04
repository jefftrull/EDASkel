// Test of PRIMA model reduction using Eigen
// test circuit is same as previous odeint+Eigen(MNA) RLC tank circuit
// Jeff Trull 2013-06-01

#include <vector>
#include <algorithm>
using namespace std;
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric;
#include <Eigen/Dense>
#include <Eigen/StdVector>
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
  Matrix3d coeff_;
  Matrix<double, 3, 1> input_;   // for applying excitation to state variables
  Matrix<double, 2, 3> output_;  // for generating output from state variables

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
    stamp_i(G, 1, 3);   // connect current derivative to associated voltage

    // Rewrite state equation as C*dx/dt = -G*x + B*u(t)
    // Plus output v(t) = L.transpose() * x

    // Now reduce C, G, B, and L so they refer to a reduced number of states,
    // but their interface (inputs and outputs) remains the same

    // PRIMA

    // Step 1: create B and L (input and output) matrices

    // We have one input: the applied voltage source, so u(t) is a 1x1 vector
    // The matrix B, when multiplied by u(t), must produce a column vector with
    // -Va in the third row (finishing out the BCE 0 = V0 - Va)
    Matrix<double, 4, 1> B; B << 0,  0, -1, 0;

    // We have two (system) outputs: the current of the voltage input, and the
    // voltage on the second (LC tank) node.  The matrix L, transposed, times
    // the system variables should produce these values.  In other words:
    // L.transpose() * x = [Ia, V1] (a column vector)
    // So L.transpose() is 2x4, and L is thus 4x2
    Matrix<double, 4, 2> L; L << 0, 0
                               , 0, -1  // for extracting V1 (BOZO seems wrong)
                               , -1, 0  // for extracting Ia
                               , 0, 0 ;

    // Step 2: Solve GR = B for R
    Matrix<double, 4, 1> R = G.colPivHouseholderQr().solve(B);

    // "X" is quite overloaded in the PRIMA paper.  I have to name them uniquely.
    // Variable translation key:
    // Just plain X is the X[k] matrices that are kept throughout the algorithm
    // Inside the loop Xk[j] stores the temporary sequence that uses parenthesized (j) superscripts
    // Xfinal is the result matrix we get at the end

    typedef Matrix<double, 4, Dynamic> Matrix4dX;
    typedef aligned_allocator<Matrix4dX> Allocator4dX;
    typedef vector<Matrix4dX, Allocator4dX> Matrix4dXList;
    Matrix4dXList X;  // each item in the sequence of Xk's may have a different column count

    // Step 3: Set X0 to the orthonormal basis of R as determined by QR factorization
    // As pointed out in Boyer for this case (R is a column vector) we simply normalize
    X.push_back(R.normalized());  // once normalized, is its own orthonormal basis

    // Step 4: Set n = floor(q/N)+1 if q/N is not an integer, and q/N otherwise
    // q = desired number of state variables in the reduced system = 3 (inductor current,
    //     capacitor voltage, input voltage - or some analogs for those)
    // N = number of ports, 2 in our case
    size_t n = 2;

    // Step 5: Block Arnoldi (see Boyer for detailed explanation)
    for (size_t k = 1; k <= n; ++k)
    {
      vector<Matrix<double, 4, 1>,
             aligned_allocator<Matrix<double, 4, 1> > > Xk(n+1);

       // set V = C * X[k-1]
       auto V = C * X[k-1];

       // solve G*X[k][0] = V for X[k][0]
       Xk[0] = G.colPivHouseholderQr().solve(V);

       for (size_t j = 1; j <= k; ++j)
       {
          // H = X[k-j].transpose() * X[k][j-1]
          auto H = X[k-j].transpose() * Xk[j-1];

          // X[k][j] = X[k][j-1] - X[k-j]*H
          Xk[j] = Xk[j-1] - X[k-j] * H;
       }

       // set X[k] to the orthonormal basis of X[k][k] via QR factorization
       X.push_back(Xk[k].normalized()) ;   // a single column, cannot factorize...


    }

    // Step 6: Set Xfinal to the concatenation of all those bases we calculated above,
    //         truncated to q==3 columns
    size_t cols = accumulate(X.begin(), X.end(), 0,
                             [](size_t sum,
                                Matrix4dX const& m) { return sum + m.cols(); });

    Matrix<double, 4, 3> Xfinal(4, cols);
    size_t col = 0;
    for (size_t k = 0; (k <= n) && (col < cols); ++k)
    {
      // copy columns from X[k] to Xfinal
      for (int j = 0; (j < X[k].cols()) && (col < cols); ++j)
      {
        Xfinal.col(col++) = X[k].col(j);
      }
    }

    // Step 7: Compute the C and G matrices in the new state variables using X
    MatrixXd Cprime = Xfinal.transpose() * C * Xfinal;
    MatrixXd Gprime = Xfinal.transpose() * G * Xfinal;
    auto     Bprime = Xfinal.transpose() * B;

    // We are done with PRIMA.  Set up suitable matrices so we can calculate
    // new state and output values:

    // this array multiplies the state vector
    // BOZO review this decomposition for correctness (and the next one too)
    coeff_ = Cprime.ldlt().solve(-1.0 * Gprime);
     
    // this array multiplies the input vector (which is a 1x1)
    input_ = Cprime.ldlt().solve(Bprime);

    // now the new state equation is dx/dt = coeff_ * x + input_ * Va;

    // Similarly we have to be able to translate our new state variables into
    // our Ia and V1 outputs
    auto Lprime = Xfinal.transpose() * L;

    // and the state-to-output translation is just Lprime
    output_ = Lprime.transpose();

  }

  void operator() (const state_type x, state_type& dxdt, double t) {

    double Va = 0;
    if (t > 0.0)
    {
      Va = vin_;   // step function at time 0
    }

    Map<const Matrix<double, 3, 1> > xvec(x.data());  // turn state vector into Eigen matrix
    Matrix<double, 1, 1> u; u << Va;                  // input excitation is 1x1 vector of Va

    Map<Matrix<double, 3, 1> > result(dxdt.data());
    result = coeff_ * xvec + input_ * u;

  }

  Matrix<double, 2, 3> output() const
  {
    return output_;
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


int main() {
  // pick some seemingly sensible values
  double v = 3.3;
  double r = 100;
  double l = 20.0e-6;   // 20 uH
  double c = 20.0e-9;   // 20 nF
  rlc_tank ckt(r, l, c, v);


  // initial state: all voltages and currents 0
  state_type x({0.0, 0.0, 0.0});

  vector<state_type> state_history;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 10e-6, 1e-6,
                     push_back_state_and_time( state_history, times ) );

  auto output_translator = ckt.output();

  for (size_t i = 0; i < times.size(); ++i) {
    // format for gnuplot.  Only vin and vout are interesting
    Map<const Matrix<double, 3, 1> > statevec(state_history[i].data());
    double vout = (output_translator * statevec)(1);

    cout << times[i] << " " << vout << endl;
  }

  return 0;

}
