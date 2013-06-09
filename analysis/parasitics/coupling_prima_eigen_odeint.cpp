// Test of odeint using intersignal coupling
// Now with the use of Eigen to separate the equations
// and to perform PRIMA model reduction
// Jeff Trull 2013-06-04

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

   matrix(vnodeno, istateno) = 1;   // current is taken *into* inductor or vsource
   matrix(istateno, vnodeno) = -1;
}

typedef vector<double> state_type;

struct signal_coupling {
  typedef Matrix<double, 10, 10> Matrix10d;   // for the MNA-based system

  Matrix<double, 4, 4> coeff_;                // reduced system state evolution
  Matrix<double, 4, 2> input_;                // inputs (agg/vic) to reduced system state
  Matrix<double, 2, 4> output_;               // reduced system state to current outputs

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

      // PRIMA
      
      // Aiming for 4 state variables, a significant reduction from 10...

      // Step 1: create B and L (input and output) matrices

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
      Matrix<double, 10, 2> L; L << 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , -1, 0    // extract V3 (aggressor rcvr)
                                  , 0, 0
                                  , 0, 0
                                  , 0, 0
                                  , 0, -1    // extract v7 (victim rcvr)
                                  , 0, 0
                                  , 0, 0 ;
                                  
      // Step 2: Solve GR = B for R
      Matrix<double, 10, 2> R = G.colPivHouseholderQr().solve(B);

      // set up types for various versions of "X" variables used in PRIMA
      // matrices we are gathering will all be 10 rows tall but some unknown number
      // of columns, depending on how many bases we harvest from each Krylov matrix
      typedef Matrix<double, 10, Dynamic> Matrix10dX;
      // we want to use these in std::vectors, which requires a special allocator
      // in order to be Eigen-compatible:
      typedef aligned_allocator<Matrix10dX> Allocator10dX;
      typedef vector<Matrix10dX, Allocator10dX> Matrix10dXList;
      Matrix10dXList X;   // one entry per value of "k", gathered at the end into a single Xfinal

      // Step 3: Set X[0] to the orthonormal basis of R as determined by QR factorization
      auto rQR = R.colPivHouseholderQr();
      Matrix10dX rQ = rQR.householderQ();    // gets 10x10 matrix, only part of which we need
      X.push_back(rQ.leftCols(rQR.rank()));  // only use the basis part

      std::cout << "X[0] is of rank " << R.colPivHouseholderQr().rank() << endl;
      std::cout << "X[0] has " << X[0].cols() << " columns\n";

      // Step 4: Set n = floor(q/N)+1 if q/N is not an integer, and q/N otherwise
      const size_t q = 4; // desired number of state variables in the reduced system (guessing here)
      const size_t N = 4; // number of ports.  But for us, two are in, two out.  Hm...
      size_t n = 2;       // this seems to work best despite not matching PRIMA equation :(

      // Step 5: Block Arnoldi (see Boyer for detailed explanation)
      for (size_t k = 1; k <= n; ++k)
      {
         // because X[] will vary in number of columns, so will Xk[]
         vector<Matrix<double, 10, Dynamic>,
                aligned_allocator<Matrix<double, 10, Dynamic> > > Xk(n+1);

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
         if (Xk[k].cols() == 1)
         {
            // a single column is automatically orthogonalized; just normalize
            X.push_back(Xk[k].normalized());
            std::cout << "produced 1 column on this iteration:\n" << X[k] << endl;
         } else {
            auto xkkQR = Xk[k].colPivHouseholderQr();
            Matrix10dX xkkQ = xkkQR.householderQ();
            X.push_back(xkkQ.leftCols(xkkQR.rank()));
            std::cout << "produced " << xkkQR.rank() << " columns on this iteration:\n" << X[k] << endl;
         }
      }

      // Step 6: Set Xfinal to the concatenation of all those bases we calculated above,
      //         truncated to q columns
      size_t cols = accumulate(X.begin(), X.end(), 0,
                               [](size_t sum, Matrix10dX const& m) { return sum + m.cols(); });
      cols = std::min(q, cols);  // truncate to q

      Matrix10dX Xfinal(10, cols);
      size_t col = 0;
      for (size_t k = 0; (k <= n) && (col < cols); ++k)
      {
         // copy columns from X[k] to Xfinal
         for (int j = 0; (j < X[k].cols()) && (col < cols); ++j)
         {
            Xfinal.col(col++) = X[k].col(j);
         }
      }

      std::cout << "Xfinal is:\n" << Xfinal << endl;
      std::cout << "Xfinal transpose times Xfinal is:\n" << Xfinal.transpose() * Xfinal << endl;

      // Step 7: Compute the C and G matrices in the new state variables using X
      MatrixXd Cprime = Xfinal.transpose() * C * Xfinal;
      MatrixXd Gprime = Xfinal.transpose() * G * Xfinal;
      auto     Bprime = Xfinal.transpose() * B;

      // PRIMA done.  Set up reduced system matrices

      coeff_ = Cprime.ldlt().solve(-1.0 * Gprime);   // state evolution
      input_ = Cprime.ldlt().solve(Bprime);          // from inputs to state variables
      auto Lprime = Xfinal.transpose() * L;
      output_ = Lprime.transpose();                  // from state to outputs

      std::cout << "Cprime:\n" << Cprime << endl;
      std::cout << "Gprime:\n" << Gprime << endl;
      std::cout << "Bprime:\n" << Bprime << endl;

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

      std::cout << "calculated real_start=" << real_start << " and real_end=" << real_end << endl;

      if (t >= real_end) {
        Vagg = (agg_slew_ < 0) ? 0 : v_;      // straight to final value
      } else if (t <= real_start) {
        Vagg = (agg_slew_ > 0) ? 0 : v_;      // remain at initial value
      } else {
        Vagg += agg_slew_ * (t - real_start);  // proportional to time in ramp
      }

      std::cout << "so for t=" << t << " Vagg=" << Vagg << endl;

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
        
    Map<const Matrix<double, 4, 1> > xvec(x.data());  // turn state vector into Eigen matrix
    Matrix<double, 2, 1> u; u << Vagg, Vvic;          // input excitation

    Map<Matrix<double, 4, 1> > result(dxdt.data());
    result = coeff_ * xvec + input_ * u;              // sets dxdt via reference

    std::cout << "at time " << t << " with a state vector of\n" << xvec << endl << " and inputs of\n";
    std::cout << u << endl << "we applied coefficients of " << coeff_ << endl;
    std::cout << " and input applicator of " << input_ << endl;
    std::cout << " to come up with a dx/dt value of " << result << endl;

  }

  Matrix<double, 2, 4> output() const
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
  state_type x({0.0, 0.0, 0.0, 0.0});

  vector<state_type> state_history;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 1000e-12, 1e-12,
                     push_back_state_and_time( state_history, times ) );

  auto output_translator = ckt.output();

  for (size_t i = 0; i < times.size(); ++i) {
    // format for gnuplot.  We are remembering the output voltages
    Map<const Matrix<double, 4, 1> > statevec(state_history[i].data());
    double vagg = (output_translator * statevec)(0);
    double vvic = (output_translator * statevec)(1);

    cout << times[i] << " " << vagg << " " << vvic << endl;

  }

  // find the highest voltage on the victim (which is supposed to be low)
  cerr << "max victim excursion is: " << max_voltage(state_history, 7) << endl;

  cerr << "driver delay is: " << delay(times, state_history, v/2.0, 1, 3, RiseRise) << endl;

  return 0;

}
