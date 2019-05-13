// Tests for Model Order Reduction algorithms
// Copyright (C) 2014 Jeffrey Elliot Trull <edaskel@att.net>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#define BOOST_TEST_MODULE Model Order Reduction
#include <boost/test/included/unit_test.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../mna.hpp"

BOOST_AUTO_TEST_CASE( RC_3x3 ) {
    // A simple resistive trace with a single resistor and capacitor
    using Eigen::Matrix;
    using namespace EDASkel::analysis::mna;

    const size_t statecnt = 3;  // Vinput, Voutput, Iinput
    typedef Matrix<double, statecnt, statecnt> eqn_matrix_t;
    eqn_matrix_t G = eqn_matrix_t::Zero();
    eqn_matrix_t C = eqn_matrix_t::Zero();
    stamp(G, 0, 1, 0.01f);      // 100ohm from voltage to load (conductance)
    stamp(C, 1, 1e-12);         // 1pF load
    stamp_i(G, 0, 2);           // connect voltage source to its output current
    
    Matrix<double, statecnt, 1> B, D;
    B << 0,  0, -1;             // voltage source is input
    D << 0,  1,  0;             // end of trace is output

    Matrix<double, Dynamic, Dynamic> Greg, Creg;
    Matrix<double, Dynamic, 1> Breg, Dreg;
    Matrix<double, 1, 1> E;
    std::tie(Greg, Creg, Breg, Dreg, E) = regularize_natarajan(G, C, B, D);

    BOOST_CHECK( !isSingular(Creg) );

    typedef MatrixVector<double, 1, 1> moments_t;
    moments_t moments_nr = moments(Greg, Creg, Breg, Dreg, E, 5);

    // A single RC is simple enough to calculate moments by hand
    // They are (-RC)^i for i=0, 1, 2...
    BOOST_CHECK_CLOSE_FRACTION(      1, moments_nr[0](0, 0), 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( -1e-10, moments_nr[1](0, 0), 0.001 );
    BOOST_CHECK_CLOSE_FRACTION(  1e-20, moments_nr[2](0, 0), 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( -1e-30, moments_nr[3](0, 0), 0.001 );
    BOOST_CHECK_CLOSE_FRACTION(  1e-40, moments_nr[4](0, 0), 0.001 );

    // Now test Su regularization for the same circuit
    Matrix<double, Dynamic, Dynamic> Greg_su, Creg_su;
    Matrix<double, Dynamic, 1> Breg_su, Dreg_su;
    std::tie(Greg_su, Creg_su, Breg_su, Dreg_su) = regularize_su(G, C, B, D);

    BOOST_CHECK( !isSingular(Creg_su) );

    Matrix<double, 1, 1> Ezero = Matrix<double, 1, 1>::Zero();
    moments_t moments_su = moments(Greg_su, Creg_su, Breg_su, Dreg_su, Ezero, 5);

    BOOST_CHECK_CLOSE_FRACTION(      1, moments_su[0](0, 0), 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( -1e-10, moments_su[1](0, 0), 0.001 );
    BOOST_CHECK_CLOSE_FRACTION(  1e-20, moments_su[2](0, 0), 0.001 );
    BOOST_CHECK_CLOSE_FRACTION( -1e-30, moments_su[3](0, 0), 0.001 );
    BOOST_CHECK_CLOSE_FRACTION(  1e-40, moments_su[4](0, 0), 0.001 );
}

// Examples from the Natarajan paper

BOOST_AUTO_TEST_CASE( Sallen_Key_Filter ) {
    using Eigen::Matrix;
    using namespace EDASkel::analysis::mna;

    const size_t statecnt = 7;
    Matrix<double, statecnt, statecnt> G, C;
    G << 0,  0,  0,  0,  0,  1,  0,
         0,  1,  0,  0, -1,  0,  0,
         0,  0,  1,  0,  0,  0,  0,
         0,  0,  0,  2, -1,  0,  0,
         0, -1,  0, -1,  2,  0,  1,
         1,  0,  0,  0,  0,  0,  0,
         0,  0,  1, -1,  0,  0,  0;

    C << 1, -1,  0,  0,  0,  0,  0,
        -1,  2, -1,  0,  0,  0,  0,
         0, -1,  1,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  -0.001,  0,  0;

    Matrix<double, statecnt, 1> B, D;
    B << 0,  0,  0,  0,  0,  1,  0;
    D << 0,  0,  0,  0,  1,  0,  0;

    Matrix<double, Dynamic, Dynamic> Greg, Creg;
    Matrix<double, Dynamic, 1> Breg, Dreg;
    Matrix<double, 1, 1> E;
    std::tie(Greg, Creg, Breg, Dreg, E) = regularize_natarajan(G, C, B, D);

    BOOST_CHECK( E.isZero() );  // no feedthrough

    MatrixVector<double, 1, 1> moments_nr = moments(Greg, Creg, Breg, Dreg, E, 5);

    // Now calculate Natarajan manual calculated results described in the paper
    // Begin after step (b) is complete (last complete intermediate step shown)
    Matrix<double, 7, 7> C0;
    C0 << 1, -1,  0,  0,  0,  0,  0,
          0,  1, -1,  0,  0,  0,  0,
          0,  0,  0,  0, -0.001, 0, 0,
          0,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  0,  0;

    // G with the same manipulations applied, followed by putting G2 into row echelon form, reversed
    Matrix<double, 7, 7> G0;
    G0 << 0,  0,  0,  0,  0,  1,  0,
          0,  1,  0,  0, -1,  1,  0,
          0,  0,  1, -1,  0,  0,  0,
          1,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  2, -1,  0,  0,
          0,  1,  1,  0, -1,  1,  0,
          0, -1,  0, -1,  2,  0,  1;

    // B with the same row manipulations applied
    Matrix<double, 7, 1> B0;
    B0 << 0,  0,  0,  1,  0,  0,  0;
    Matrix<double, 1, 7> D0;             // no change
    D0 << 0,  0,  0,  0,  1,  0,  0;

    // swap columns 0 and 3 so G22 will be non-singular (puts a 1 in the 3,3 position)
    Eigen::PermutationMatrix<7, 7, std::size_t> g22_column_swap;
    g22_column_swap.setIdentity(7);
    g22_column_swap.applyTranspositionOnTheRight(0, 3);
    C0 = C0 * g22_column_swap;
    G0 = G0 * g22_column_swap;
    D0 = D0 * g22_column_swap;

    // Calculate submatrices
    Matrix<double, 3, 3> G11 = G0.topLeftCorner(3, 3);
    Matrix<double, 4, 3> G21 = G0.bottomLeftCorner(4, 3);
    Matrix<double, 3, 4> G12 = G0.topRightCorner(3, 4);
    Matrix<double, 4, 4> G22 = G0.bottomRightCorner(4, 4);

    Matrix<double, 3, 3> C11 = C0.topLeftCorner(3, 3);
    Matrix<double, 3, 4> C12 = C0.topRightCorner(3, 4);

    Matrix<double, 3, 1> B01 = B0.topRows(3);
    Matrix<double, 4, 1> B02 = B0.bottomRows(4);

    Matrix<double, 1, 3> D01 = D0.leftCols(3);
    Matrix<double, 1, 4> D02 = D0.rightCols(4);

    assert(!isSingular(G22));
    auto G22_LU = G22.fullPivLu();

    Matrix<double, 3, 3> G1 = G11 - G12 * G22_LU.solve(G21);
    Matrix<double, 3, 3> C1 = C11 - C12 * G22_LU.solve(G21);
    Matrix<double, 3, 1> B1 = B01 - G12 * G22_LU.solve(B02);
    Matrix<double, 3, 1> B2 =     - C12 * G22_LU.solve(B02);

    // Do Chen variable substitution
    Matrix<double, 3, 1>       Br = B1  - G1 * C1.inverse() * B2;
    Matrix<double, 3, 1>       D1 = (D01 - D02 * G22_LU.solve(G21)).transpose();
    Matrix<double, 1, 1>    Ezero = Matrix<double, 1, 1>::Zero();
    MatrixVector<double, 1, 1> moments_manual = moments(G1, C1, Br, D1, Ezero, 5);

    // Perform comparisons
    BOOST_CHECK_SMALL( moments_nr[0](0, 0), 1e-10 );   // these first two are theoretically 0
    BOOST_CHECK_SMALL( moments_nr[1](0, 0), 1e-10 );
    BOOST_CHECK_CLOSE_FRACTION( moments_manual[2](0, 0), moments_nr[2](0, 0), 1e-10 );
    BOOST_CHECK_CLOSE_FRACTION( moments_manual[3](0, 0), moments_nr[3](0, 0), 1e-10 );
    BOOST_CHECK_CLOSE_FRACTION( moments_manual[4](0, 0), moments_nr[4](0, 0), 1e-10 );

}

BOOST_AUTO_TEST_CASE( Passive_Network ) {
    using Eigen::Matrix;
    using namespace EDASkel::analysis::mna;

    const size_t statecnt = 7;
    Matrix<double, statecnt, statecnt> G, C;
    G << 0.5, -0.5,  0,    0,    0,   1,  0,
        -0.5,  0.5,  0,    0,    0,   0,  0,
         0,    0,    0.25, 0,    0,   0,  0,
         0,    0,    0,    0.2, -0.2, 0,  0,
         0,    0,    0,   -0.2,  0.2, 0,  1,
         1,    0,    0,    0,    0,   0,  0,
         0,    0,    0,    0,    1,   0,  0;

    C << 0,  0,  0,  0,  0,  0,  0,
         0,  3, -2, -1,  0,  0,  0,
         0, -2,  6, -4,  0,  0,  0,
         0, -1, -4,  5,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0, -5;

    Matrix<double, statecnt, 1> B;
    B << 0,  0,  0,  0,  0,  1,  0;
    Matrix<double, statecnt, 2> D;
    D << 0,  0,
         0,  1,
         0,  0,
         0, -1,
         1,  0,
         0,  0,
         0,  0;

    Matrix<double, Dynamic, Dynamic> Greg, Creg;
    Matrix<double, Dynamic, 1> Breg;
    Matrix<double, Dynamic, 2> Lreg;
    Matrix<double, 2, 1> E;
    std::tie(Greg, Creg, Breg, Lreg, E) = regularize_natarajan(G, C, B, D);
    MatrixVector<double, 2, 1> moments_nr = moments(Greg, Creg, Breg, Lreg, E, 2);

    // values from the paper
    Matrix<double, 3, 3> Aexpected;   // -C^-1 * G
    Aexpected <<
        -0.267, -0.500, -0.292,
        -0.267, -0.714, -0.327,
        -0.267, -0.571, -0.345 ;      // <<< error in paper! not .245...
    Matrix<double, 3, 1> Bexpected;   //  C^-1 * B1
    Bexpected << 0.5, 0.714, 0.571;
    Matrix<double, 3, 1> B1expected;  //  C^-1 * B2
    B1expected << 0.667, 0.667, 0.667;
    Matrix<double, 2, 3> Dexpected;
    Dexpected <<
         1, 2.5, 1.25,
        -1, 1,   0;
    Matrix<double, 2, 1> Eexpected;
    Eexpected << -2.5, 0;

    // These results are of the form dX/dt = A*X + B*u + B1*du/dt
    // To get moments, use Y = D*X + E, then take Laplace transform
    // I get H(s) = D * (B + sB1) / (sI - A) + E
    // H(0)  = -D * A^-1 * B + E
    Matrix<double, Dynamic, 1> Rexpected = Aexpected.fullPivHouseholderQr().solve(Bexpected);
    Matrix<double, 2, 1> moment_expected0 = -Dexpected * Rexpected + Eexpected;
    BOOST_CHECK_SMALL( moment_expected0(0, 0), 1e-10 );

    // H'(0) = -D * A^-2 * (B + A * B1)
    Matrix<double, 2, 1> moment_expected1 = -Dexpected *
        Aexpected.fullPivHouseholderQr().solve(
            Aexpected.fullPivHouseholderQr().solve(Bexpected + Aexpected * B1expected));

    // and if we needed to do additional ones:
    // H''(0) = 2 * D * A^-3 * (B + A * B1)
    // M2 = H''(0) / (2!)
    // etc.

    // Check vs. values from paper.  Using a fairly loose tolerance b/c paper only gives 3 digits of precision
    BOOST_CHECK_SMALL( moments_nr[0](0, 0), 1e-10 );
    BOOST_CHECK_CLOSE_FRACTION( moment_expected0(1, 0), moments_nr[0](1, 0), 0.01 );
    BOOST_CHECK_SMALL( moments_nr[1](0, 0), 1e-10 );
    BOOST_CHECK_CLOSE_FRACTION( moment_expected1(1, 0), moments_nr[1](1, 0), 0.01 );

    // do the same for Su
    Matrix<double, Dynamic, Dynamic> Greg_su, Creg_su;
    Matrix<double, Dynamic, 1> Breg_su;
    Matrix<double, Dynamic, 2> Lreg_su;
    Matrix<double, 2, 1> Ezero = Matrix<double, 2, 1>::Zero();
    std::tie(Greg_su, Creg_su, Breg_su, Lreg_su) = regularize_su(G, C, B, D);
    MatrixVector<double, 2, 1> moments_su = moments(Greg_su, Creg_su, Breg_su, Lreg_su, Ezero, 2);

    BOOST_CHECK_SMALL( moments_su[0](0, 0), 1e-10 );
    BOOST_CHECK_CLOSE_FRACTION( moment_expected0(1, 0), moments_su[0](1, 0), 0.01 );
    BOOST_CHECK_SMALL( moments_su[1](0, 0), 1e-10 );
    BOOST_CHECK_CLOSE_FRACTION( moment_expected1(1, 0), moments_su[1](1, 0), 0.01 );

}
