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
  static const size_t q = 6; // desired number of state variables. 10 is the "natural" count.

  MatrixXd coeff_;           // reduced system state evolution
  MatrixXd input_;           // inputs (agg/vic) to reduced system state
  MatrixXd output_;          // reduced system state to chosen outputs

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
      // aggressor and victim driver source currents will be nodes 8 and 9
      stamp_i(G, 0, 12);
      stamp_i(G, 6, 13);

      // PRIMA
      // This is a famous model reduction technique invented around 1997/8 at CMU
      // I am generally following the treatment in Odabasioglu, IEEE TCAD, August 1998
      // They base their work on the prior Block Arnoldi algorithm, for a helpful
      // explanation of which you can find in their 6th citation:
      // D. L. Boley, “Krylov space methods on state-space control models”
      
      // Aiming for q state variables, a significant reduction from 10...

      // Step 1: create B and L (input and output) matrices

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

      // Four outputs, for viewing and performing measurements
      Matrix<double, 14, 4> L; L << 0, 0, 0, 0
                                  , 1, 0, 0, 0    // extract V1 (aggressor driver output)
                                  , 0, 0, 0, 0
                                  , 0, 0, 0, 0
                                  , 0, 0, 0, 0
                                  , 0, 0, 0, 1    // extract V5 (aggressor receiver)
                                  , 0, 0, 0, 0
                                  , 0, 0, 0, 0
                                  , 0, 1, 0, 0    // extract v8 (victim coupling node)
                                  , 0, 0, 0, 0
                                  , 0, 0, 0, 0
                                  , 0, 0, 1, 0    // extract v11 (victim rcvr)
                                  , 0, 0, 0, 0
                                  , 0, 0, 0, 0 ;
                                  
      // Step 2: Solve GR = B for R
      Matrix<double, 14, 2> R = G.fullPivHouseholderQr().solve(B);

      // set up types for various versions of "X" variables used in PRIMA
      // matrices we are gathering will all be 14 rows tall but some unknown number
      // of columns, depending on how many bases we harvest from each Krylov matrix
      typedef Matrix<double, 14, Dynamic> Matrix14dX;
      // we want to use these in std::vectors, which requires a special allocator
      // in order to be Eigen-compatible:
      typedef aligned_allocator<Matrix14dX> Allocator14dX;
      typedef vector<Matrix14dX, Allocator14dX> Matrix14dXList;
      Matrix14dXList X;   // one entry per value of "k", gathered at the end into a single Xfinal

      // Step 3: Set X[0] to the orthonormal basis of R as determined by QR factorization
      auto rQR = R.fullPivHouseholderQr();
      Matrix14dX rQ = rQR.matrixQ();    // gets 14x14 matrix, only part of which we need
      X.push_back(rQ.leftCols(rQR.rank()));  // only use the basis part

      // Step 4: Set n = floor(q/N)+1 if q/N is not an integer, and q/N otherwise
      const size_t N = 4; // number of ports.  But for us, two are in, two out.  Hm...
      size_t n = (q % N) ? (q/N + 1) : (q/N);

      // Step 5: Block Arnoldi (see Boyer for detailed explanation)
      for (size_t k = 1; k <= n; ++k)
      {
         // because X[] will vary in number of columns, so will Xk[]
         vector<Matrix<double, 14, Dynamic>,
                aligned_allocator<Matrix<double, 14, Dynamic> > > Xk(n+1);  // note to self: the size seems wrong (too large)

         // set V = C * X[k-1]
         auto V = C * X[k-1];

         // solve G*X[k][0] = V for X[k][0]
         Xk[0] = G.fullPivHouseholderQr().solve(V);

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
         } else {
            auto xkkQR = Xk[k].fullPivHouseholderQr();
            Matrix14dX xkkQ = xkkQR.matrixQ();
            X.push_back(xkkQ.leftCols(xkkQR.rank()));
         }
      }

      // Step 6: Set Xfinal to the concatenation of all those bases we calculated above,
      //         truncated to q columns
      size_t cols = accumulate(X.begin(), X.end(), 0,
                               [](size_t sum, Matrix14dX const& m) { return sum + m.cols(); });
      cols = std::min(q, cols);  // truncate to q

      Matrix14dX Xfinal(14, cols);
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

      // PRIMA done.  Next step: Regularize inputs (see RLC example for more detail)
      coeff_ = Cprime.ldlt().solve(-1.0 * Gprime);   // state evolution
      input_ = Cprime.ldlt().solve(Bprime);          // from inputs to state variables
      Matrix<double, Dynamic, 4> Lprime = Xfinal.transpose() * L;
      output_ = Lprime.transpose();                  // from state to outputs

      // Compare block moments between original and reduced model
      // Prima claims to produce the same moments up to floor(q/N)
      Matrix<double, 14, 14> A = C.fullPivHouseholderQr().solve(-1.0 * G);
      Matrix<double, q, q> Aprime = coeff_;
      Matrix<double, q, 2> Rprime = Gprime.fullPivHouseholderQr().solve(Bprime);
      Matrix<double, 14, 14> AtotheI = Matrix<double, 14, 14>::Identity();
      Matrix<double, q, q> AprimetotheI = Matrix<double, q, q>::Identity();
      for (size_t i = 0; i < (q/N); ++i) {
        std::cerr << "moment " << i << " of original model is:" << std::endl;
        std::cerr << L.transpose() * AtotheI * R << std::endl;
        std::cerr << "moment " << i << " of reduced model is:" << std::endl;
        std::cerr << Lprime.transpose() * AprimetotheI * Rprime << std::endl;
        AtotheI = A * AtotheI;
        AprimetotheI = Aprime * AprimetotheI;
      }
      
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
        
    Map<const Matrix<double, q, 1> > xvec(x.data());  // turn state vector into Eigen matrix
    Matrix<double, 2, 1> u; u << Vagg, Vvic;          // input excitation

    Map<Matrix<double, q, 1> > result(dxdt.data());
    result = coeff_ * xvec + input_ * u;              // sets dxdt via reference

  }

  Matrix<double, 4, q> output() const
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
  state_type x(signal_coupling::q, 0.0);

  vector<state_type> state_history;
  vector<state_type> outputs;
  vector<double>     times;

  odeint::integrate( ckt, x, 0.0, 1000e-12, 1e-12,
                     push_back_state_and_time( state_history, times ) );

  // generate observable (output) values by applying output transform matrix to each state
  auto const output_translator = ckt.output();
  transform(state_history.begin(), state_history.end(), back_inserter(outputs),
            [&output_translator]
            (const state_type& state) -> state_type {
              // access vector state data through an Eigen map so we can multiply
              Map<const Matrix<double, signal_coupling::q, 1> > statevec(state.data());
              state_type ovec(4);   // result
              Map<Matrix<double, 4, 1> > outputvec(ovec.data(), 4, 1);   // also needs a Map
              outputvec = output_translator * statevec;  // multiply into vector via Map
              return ovec;
            });

  for (size_t i = 0; i < times.size(); ++i) {
    // format for gnuplot.  We are remembering the output voltages
    // time, vagg, vcoup, vvic (vagg_rcv omitted)
    cout << times[i] << " " << outputs[i][0] << " " << outputs[i][1] << " " << outputs[i][2] << endl;
  }

  // find the highest voltage on the victim (which is supposed to be low)
  cerr << "max victim excursion is: " << max_voltage(outputs, 2) << endl;

  cerr << "driver delay is: " << delay(times, outputs, v/2.0, 0, 3, RiseRise) << endl;

  return 0;

}

const size_t signal_coupling::q;
