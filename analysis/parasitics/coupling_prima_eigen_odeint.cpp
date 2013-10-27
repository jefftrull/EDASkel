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
#include <Eigen/Eigenvalues>
using namespace Eigen;
 
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
  static const size_t q = 8; // desired number of state variables in the reduced system
                             // 12 is the "natural" count (1 per circuit node)
  static const size_t N = 3; // port count (one per input/output: two drivers and the victim rcvr)

  // final system for simulation: 
  Matrix<double, q+1, q+1> coeff_;        // reduced system state evolution
  Matrix<double, q+1, 2> input_;          // inputs (agg/vic) to reduced system state
  Matrix<double, 1, q+1> output_;         // reduced system state to chosen outputs

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
      typedef Matrix<double, 15, 15> Matrix15d;
      Matrix15d C = Matrix15d::Zero(), G = Matrix15d::Zero();

      stamp(G, 0, 1, 1/agg_imp);   // driver impedances
      stamp(G, 6, 7, 1/vic_imp);

      stamp(C, 1, agg_c1 / 4.0);   // beginning of first stage pi model
      stamp(C, 7, vic_c1 / 4.0);
      stamp(G, 1, 2, 2/agg_r1);    // first stage pi resistance
      stamp(G, 7, 8, 2/vic_r1);
      stamp(C, 2, agg_c1 / 4.0);   // end of first stage pi model
      stamp(C, 8, vic_c1 / 4.0);

      stamp(C, 2, agg_c1 / 4.0);   // beginning of second stage pi model
      stamp(C, 8, vic_c1 / 4.0);
      stamp(G, 2, 3, 2/agg_r1);    // second stage pi resistance
      stamp(G, 8, 9, 2/vic_r1);
      stamp(C, 3, agg_c1 / 4.0);   // end of second stage pi model
      stamp(C, 9, vic_c1 / 4.0);

      stamp(C, 3, 9, coup_c);      // coupling capacitance between traces

      stamp(C, 3, agg_c2 / 4.0);   // beginning of third stage pi model
      stamp(C, 9, vic_c2 / 4.0);
      stamp(G, 3, 4, 2/agg_r2);    // third stage pi resistance
      stamp(G, 9, 10, 2/vic_r2);
      stamp(C, 4, agg_c2 / 4.0);   // end of third stage pi model
      stamp(C, 10, vic_c2 / 4.0);

      stamp(C, 4, agg_c2 / 4.0);   // beginning of fourth stage pi model
      stamp(C, 10, vic_c2 / 4.0);
      stamp(G, 4, 5, 2/agg_r2);    // fourth stage pi resistance
      stamp(G, 10, 11, 2/vic_r2);
      stamp(C, 5, agg_c2 / 4.0);   // end of fourth stage pi model
      stamp(C, 11, vic_c2 / 4.0);

      stamp(C, 5, agg_cl);
      stamp(C, 11, vic_cl);

      // add an additional pair of equations for the independent sources
      // aggressor and victim driver source currents will be nodes 12 and 13
      stamp_i(G, 0, 12);
      stamp_i(G, 6, 13);
      stamp_i(G, 11, 14);           // victim receiver "driver" (we won't use it) current

      // PRIMA
      // This is a famous model reduction technique invented around 1997/8 at CMU
      // I am generally following the treatment in Odabasioglu, IEEE TCAD, August 1998
      // They base their work on the prior Block Arnoldi algorithm, a helpful
      // explanation of which you can find in their 6th citation:
      // D. L. Boley, “Krylov space methods on state-space control models”
      
      // Aiming for q state variables, a significant reduction from 15...

      // Step 1: create B and L (input and output) matrices

      // In the interest of following the PRIMA paper as closely as possible, we will
      // have three generic "ports" that are theoretically attached to voltage sources.
      // In fact two of these will be driven through drivers, and the last will be attached
      // to a load (the victim receiver).  But for PRIMA, they are all I/O (voltage in, current out)

      // connecting inputs
      Matrix<double, 15, N> B; B << 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , 0, 0, 0
                                  , -1, 0, 0   // insert Vagg = V0
                                  , 0, -1, 0   // insert Vvic = V6
                                  , 0, 0, -1 ; // insert Vvicrcv = V11

      // And the same nodes are outputs:
      Matrix<double, 15, N> L = -B;
                                  
      // Step 2: Solve GR = B for R
      Matrix<double, 15, N> R = G.fullPivHouseholderQr().solve(B);

      // set up types for various versions of "X" variables used in PRIMA
      // matrices we are gathering will all be 15 rows tall but some unknown number
      // of columns, depending on how many bases we harvest from each Krylov matrix

      // we want to use these in std::vectors, which requires a special allocator
      // in order to be Eigen-compatible:
      typedef aligned_allocator<Matrix15dX> Allocator15dX;
      typedef vector<Matrix15dX, Allocator15dX> Matrix15dXList;
      Matrix15dXList X;   // one entry per value of "k", gathered at the end into a single Xfinal

      // Step 3: Set X[0] to the orthonormal basis of R as determined by QR factorization
      auto rQR = R.fullPivHouseholderQr();
      Matrix15dX rQ = rQR.matrixQ();    // gets 15x15 matrix, only part of which we need
      X.push_back(rQ.leftCols(rQR.rank()));  // only use the basis part

      // Step 4: Set n = floor(q/N)+1 if q/N is not an integer, and q/N otherwise
      size_t n = (q % N) ? (q/N + 1) : (q/N);

      // Step 5: Block Arnoldi (see Boyer for detailed explanation)
      for (size_t k = 1; k <= n; ++k)
      {
         // because X[] will vary in number of columns, so will Xk[]
         vector<Matrix<double, 15, Dynamic>,
                aligned_allocator<Matrix<double, 15, Dynamic> > > Xk(k+1);

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
            Matrix15dX xkkQ = xkkQR.matrixQ();
            X.push_back(xkkQ.leftCols(xkkQR.rank()));
         }


      }

      // Step 6: Set Xfinal to the concatenation of X[0] to X[n-1],
      //         truncated to q columns
      // Note this implies we calculated X[n] unnecessarily... I'm not sure what to make of this
      size_t cols = accumulate(X.begin(), X.end()-1, 0,
                               [](size_t sum, Matrix15dX const& m) { return sum + m.cols(); });
      cols = std::min(q, cols);  // truncate to q

      Xfinal = Matrix15dX(15, cols);
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
      Matrix<double, q, q> Cprime = Xfinal.transpose() * C * Xfinal;  // TODO: look at SVD here
      Matrix<double, q, q> Gprime = Xfinal.transpose() * G * Xfinal;
      Matrix<double, q, N> Bprime = Xfinal.transpose() * B;
      Matrix<double, q, N> Lprime = Xfinal.transpose() * L;

      // evaluate results

      // Compare block moments between original and reduced model
      // Prima claims to produce the same moments up to floor(q/N)
      Matrix<double, 15, 15> A = G.fullPivHouseholderQr().solve(C);
      auto GprimeQR = Gprime.fullPivHouseholderQr();
      Matrix<double, q, q> Aprime = GprimeQR.solve(Cprime);
      Matrix<double, q, N> Rprime = GprimeQR.solve(Bprime);
      Matrix<double, 15, 15> AtotheI = Matrix<double, 15, 15>::Identity();

      Matrix<double, q, q> AprimetotheI = Matrix<double, q, q>::Identity();
      for (size_t i = 0; i < (q/N); ++i) {
        std::cerr << "moment " << i << " of original model is:" << std::endl;
        std::cerr << L.transpose() * AtotheI * R << std::endl;
        std::cerr << "moment " << i << " of reduced model is:" << std::endl;
        std::cerr << Lprime.transpose() * AprimetotheI * Rprime << std::endl;
        AtotheI = A * AtotheI;
        AprimetotheI = Aprime * AprimetotheI;
      }
      
      // Also compare eigenvalues, which are evidently the reciprocals of the poles
      // They are not guaranteed to be the same but should be "similar" (?)
      std::cerr << "eigenvalues of original model are:\n" << EigenSolver<MatrixXd>(A).eigenvalues() << std::endl;
      std::cerr << "eigenvalues of reduced model are:\n" << EigenSolver<decltype(Aprime)>(Aprime).eigenvalues() << std::endl;

      // PRIMA done.  Next step: Construct "Direct Stamp" reduced realization
      // The system equation (see "(50)" in the paper) hooks up the ports of
      // the reduced network to the remainder of the system for simulation
      // Fortunately there is no remainder for us as we have presented
      // the entire thing to PRIMA.  However, we will still need to turn the
      // victim receiver load (treated by PRIMA like all ports: as a voltage source)
      // into a plain circuit node so nothing flows out of the system and we
      // can simply observe the output.

      // state variables in our direct system are:
      // 0, 1, 2: node voltages for the two drivers and the victim receiver
      //          (two are directly driven while one will "float")
      // 3, 4: source currents for the two drivers
      // 5: "source current" for the last "driver" - always 0 as this is a receiver
      // 6.. : reduced system internal state variables (of quantity q)

      // equations for our direct system will be:
      // 2 equations to connect nodes 0 and 1 to the driver voltages
      // 1 equation connecting the victim receiver current to a small load
      // 3 equations connecting the reduced state to the current "outputs"
      // q equations for the normal reduced state evolution (will use state vars 0-2
      // for the input excitation)
      typedef Matrix<double, 6+q, 6+q> direct_eqn_matrix_t;
      direct_eqn_matrix_t Gdirect = direct_eqn_matrix_t::Zero();   // multiplies state vector
      direct_eqn_matrix_t Cdirect = direct_eqn_matrix_t::Zero();   // multiplies 1st derivative of state
      Matrix<double, 6+q, 2> Bdirect;
      Bdirect.block<2, 2>(0, 0) = -Matrix<double, 2, 2>::Identity();  // connects sources to nodes
      Matrix<double, 6+q, 1> Ldirect = Matrix<double, 6+q, 1>::Zero();
      Ldirect(2, 0) = 1;              // prima port 2 voltage is the sole output

      // Generate equations of the form Cdirect*(d(state)/dt) = -Gdirect*state + Bdirect*input;

      // hook up independent sources in first two rows: connects inputs to Vagg/Vvic via Bdirect
      Gdirect.block<2, 2>(0, 0) = -Matrix<double, 2, 2>::Identity(); 

      // Fix an issue where the generation of the reduced model for simulation encounters
      // a singularity in the G22 matrix.  Probably a hack.
      Gdirect(2, 5) = 1; Cdirect(2, 2) = 1e-15;  // exposed cap at end of wire, formerly driven by voltage source

      // connect reduced state to currents (line 3 of "(50)")
      Gdirect.block<3, 3>(3, 3) = Matrix<double, 3, 3>::Identity();  // extract source currents
      Gdirect.block<3, q>(3, 6) = -Lprime.transpose();

      // connect reduced system  (line 4 of "(50)")

      auto Gprime_inv_x_Bprime = GprimeQR.solve(Bprime);
      auto Gprime_inv_x_Cprime = GprimeQR.solve(Cprime);

      Gdirect.block<q, 3>(6, 0) = Gprime_inv_x_Bprime;
      Gdirect.block<q, q>(6, 6) = Matrix<double, q, q>::Identity();
      Cdirect.block<q, q>(6, 6) = Gprime_inv_x_Cprime;
      
      // Direct Stamp realization constructed

      // now use the techniques described in Su (Proc 15th ASP-DAC, 2002) to reduce
      // this set of equations so the state variable derivatives have coefficients
      // Otherwise we cannot integrate to get the time domain result...

      // Use Eigen reductions to find zero rows
      auto zero_rows = (Cdirect.array() == 0.0).rowwise().all();   // per row "all zeros"
      std::size_t zero_count = zero_rows.count();
      std::size_t nonzero_count = 6+q - zero_count;

      direct_eqn_matrix_t permut = direct_eqn_matrix_t::Identity();   // null permutation to start
      std::size_t i, j;
      for (i = 0, j=(6+q-1); i < j;) {
        // loop invariant: rows > j are all zero; rows < i are not
        while ((i < 6+q) && !zero_rows(i)) ++i;
        while ((j > 0) && zero_rows(j)) --j;
        if (i < j) {
          // exchange rows i and j via the permutation vector
          permut(i, i) = 0; permut(j, j) = 0;
          permut(i, j) = 1; permut(j, i) = 1;
          ++i; --j;
        }
      }

      // 2. Apply permutation to MNA matrices
      direct_eqn_matrix_t CdirectP = permut * Cdirect * permut;       // permute rows and columns
      direct_eqn_matrix_t GdirectP = permut * Gdirect * permut;
      Matrix<double, 6+q, 2> BdirectP = permut * Bdirect;    // permute only rows
      Matrix<double, 6+q, 1> LdirectP = permut * Ldirect;
      
      // 3. Produce reduced equations following Su (Proc. 15th ASP-DAC, 2002)

      auto G11 = GdirectP.topLeftCorner(nonzero_count, nonzero_count);
      auto G12 = GdirectP.topRightCorner(nonzero_count, zero_count);
      MatrixXd G21 = GdirectP.bottomLeftCorner(zero_count, nonzero_count);
      MatrixXd G22 = GdirectP.bottomRightCorner(zero_count, zero_count);

      auto L1 = LdirectP.topRows(nonzero_count);
      auto L2 = LdirectP.bottomRows(zero_count);

      auto B1 = BdirectP.topRows(nonzero_count);
      auto B2 = BdirectP.bottomRows(zero_count);

      MatrixXd Cred = CdirectP.topLeftCorner(nonzero_count, nonzero_count);
      auto G22QR = G22.fullPivLu();
      MatrixXd G22invG21 = G22QR.solve(G21);
      auto G22invB2 = G22QR.solve(B2);
      auto Gred = G11 - G12 * G22invG21;

      output_ = (L1.transpose() - L2.transpose() * G22invG21);  // simplify?
      auto Bred = B1 - G12 * G22invB2;
      // assuming no "D" (direct input to output) transformation needed - we can calculate it if required

      // 4. Solve to produce a pair of matrices we can use to calculation dX/dt
      auto CredQR = Cred.fullPivHouseholderQr();
      coeff_ = CredQR.solve(-1.0 * Gred);  // state evolution
      input_ = CredQR.solve(Bred);         // input -> state
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
        
    Map<const Matrix<double, q+1, 1> > xvec(x.data());   // turn state vector into Eigen matrix

    // Input excitation:  the independent variables calculated above
    Matrix<double, 2, 1> u; u << Vagg, Vvic;

    Map<Matrix<double, q+1, 1> > result(dxdt.data());
    result = coeff_ * xvec + input_ * u;              // sets dxdt via reference

  }

  Matrix<double, 1, q+1> output() const
  {
     // supply a matrix that can be used to extract the single output from the current state
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
  state_type x(signal_coupling::q+1, 0.0);

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
               Map<const Matrix<double, signal_coupling::q+1, 1> > statevec(state.data());
              state_type ovec(1);   // result
              Map<Matrix<double, 1, 1> > outputvec(ovec.data(), 1, 1);   // also needs a Map
              outputvec = output_translator * statevec;  // multiply into vector via Map
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
