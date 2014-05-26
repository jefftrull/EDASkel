// Utility functions for circuit analysis via Modified Nodal Analysis
// Copyright (C) 2013 Jeffrey Elliot Trull <edaskel@att.net>
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

#ifndef EDASKEL_ANALYSIS_MNA_HPP
#define EDASKEL_ANALYSIS_MNA_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace EDASkel { namespace analysis { namespace mna {

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

// Overloads for when the goal is to create an Eigen SparseMatrix
// In this case we build up lists of (row, col, value) triplets to be set all at once
template<typename Float>
void stamp(typename std::vector<Eigen::Triplet<Float> >& tlist,
           std::size_t i, std::size_t j, Float g)
{
   // General stamp: conductance at [i,i] and [j,j],
   // -conductance at [i,j] and [j,i]
   // Eigen takes care of summing these for us
   tlist.emplace_back(i, i, g);
   tlist.emplace_back(j, j, g);
   tlist.emplace_back(i, j, -g);
   tlist.emplace_back(j, i, -g);

}

template<typename Float>
void stamp(typename std::vector<Eigen::Triplet<Float> >& tlist,
           std::size_t i, Float g)
{
   tlist.emplace_back(i, i, g);
}

template<typename Float>
void stamp_i(typename std::vector<Eigen::Triplet<Float> >& tlist,
             std::size_t vnodeno, std::size_t istateno)
{
   tlist.emplace_back(vnodeno, istateno, 1);
   tlist.emplace_back(istateno, vnodeno, -1);
}


template<class M>
bool isSingular(const M& m) {
  // A singular matrix has at least one zero eigenvalue
   // Use the magic of Eigen reductions to implement:
   auto eigenvalues = EigenSolver<M>(m).eigenvalues();
   return ((eigenvalues.array().real() == 0.0) &&
           (eigenvalues.array().imag() == 0.0)).any();
}

// take a circuit's linear system description in G, C, B, L form and compress it so
// the resulting C array is non-singular.  Operation depends on runtime data, so
// output array dimensions are "Dynamic"
template<int icount, int ocount, int scount, typename Float = double>
std::tuple<Matrix<Float, Dynamic, Dynamic>,   // G result
           Matrix<Float, Dynamic, Dynamic>,   // C result
           Matrix<Float, Dynamic, icount>,    // B result
           Matrix<Float, Dynamic, ocount> >   // L result
regularize_su(Matrix<Float, scount, scount> const & G,
           Matrix<Float, scount, scount> const & C,
           Matrix<Float, scount, icount> const & B,
           Matrix<Float, scount, ocount> const & L) {

   // Use the techniques described in Su (Proc 15th ASP-DAC, 2002) to reduce
   // this set of equations so the state variable derivatives have coefficients
   // Otherwise we cannot integrate to get the time domain result...

   // Use Eigen reductions to find zero rows
   auto zero_rows = (C.array() == 0.0).rowwise().all();   // per row "all zeros"
   std::size_t zero_count = zero_rows.count();
   std::size_t nonzero_count = scount - zero_count;

   // 1. Generate permutation matrix to move zero rows to the bottom
   PermutationMatrix<scount, scount, std::size_t> permut;
   permut.setIdentity();      // start with null permutation
   std::size_t i, j;
   for (i = 0, j=(scount-1); i < j;) {
      // loop invariant: rows > j are all zero; rows < i are not
      while ((i < scount) && !zero_rows(i)) ++i;
      while ((j > 0) && zero_rows(j)) --j;
      if (i < j) {
         // exchange rows i and j via the permutation vector
         permut.applyTranspositionOnTheRight(i, j);
         ++i; --j;
      }
   }

   // 2. Apply permutation to MNA matrices
   typedef Matrix<Float, scount, scount> eqn_matrix_t;
   eqn_matrix_t CP = permut * C * permut;          // permute rows and columns
   eqn_matrix_t GP = permut * G * permut;
   Matrix<Float, Dynamic, icount> BP = permut * B; // permute only rows
   Matrix<Float, Dynamic, ocount> LP = permut * L;
   
   // 3. Produce reduced equations following Su (Proc. 15th ASP-DAC, 2002)

   typedef Matrix<Float, Dynamic, Dynamic> MatrixD;
   auto G11 = GP.topLeftCorner(nonzero_count, nonzero_count);
   auto G12 = GP.topRightCorner(nonzero_count, zero_count);
   MatrixD G21 = GP.bottomLeftCorner(zero_count, nonzero_count);
   MatrixD G22 = GP.bottomRightCorner(zero_count, zero_count);

   auto L1 = LP.topRows(nonzero_count);
   auto L2 = LP.bottomRows(zero_count);

   auto B1 = BP.topRows(nonzero_count);
   auto B2 = BP.bottomRows(zero_count);

   MatrixD Cred = CP.topLeftCorner(nonzero_count, nonzero_count);

   assert(!isSingular(G22));
   auto G22QR = G22.fullPivLu();

   MatrixD G22invG21 = G22QR.solve(G21);
   auto G22invB2 = G22QR.solve(B2);
   MatrixD Gred = G11 - G12 * G22invG21;

   Matrix<Float, Dynamic, ocount> Lred = (L1.transpose() - L2.transpose() * G22invG21).transpose();
   Matrix<Float, Dynamic, icount> Bred = B1 - G12 * G22invB2;

   return std::make_tuple(Gred, Cred, Bred, Lred);
}

template<int icount, int ocount, int scount, typename Float = double>
std::tuple<Matrix<Float, Dynamic, Dynamic>,   // G result
           Matrix<Float, Dynamic, Dynamic>,   // C result
           Matrix<Float, Dynamic, icount>,    // B result
           Matrix<Float, Dynamic, ocount> >   // D result
regularize_natarajan(Matrix<Float, scount, scount> const & G,
           Matrix<Float, scount, scount> const & C,
           Matrix<Float, scount, icount> const & B,
           Matrix<Float, scount, ocount> const & D) {

    // Implements the algorithm in [Natarajan]
    // Circuits, Devices and Systems, IEE Proceedings G, June 1991

    typedef Matrix<Float, Dynamic, Dynamic> MatrixD;

    // Step 1: put C into "Row Echelon" form by performing LU factorization
    auto lu = C.fullPivLu();
    MatrixD U = lu.matrixLU().template triangularView<Upper>();
    MatrixD L = lu.matrixLU().template triangularView<UnitLower>();

    // Step 2: "perform the same elementary operations on G and B"
    // given that C = P.inverse() * L * U * Q.inverse()
    // (from source it seems that permutationP/Q is inverse)
    // then to get the new G we reverse those operations:
    auto Cprime = U;
    auto P = lu.permutationP();
    auto Q = lu.permutationQ();
    auto Gprime = L.inverse() * P * G * Q;              // rows and columns
    auto Bprime = L.inverse() * P * B;                  // rows only
    // The D input is like L in PRIMA but this algorithm uses the transpose
    auto Dprime = D.transpose() * Q;                    // columns only

    // calculate submatrices
    auto k = lu.rank();

    // Step 3: "Convert [G21 G22] matrix into row echelon form starting from the last row"
    // Since the order of the rows is irrelevant, I'll perform the decomposition, then
    // combine reversing the rows with the row reordering the LU decomposition produces

    MatrixD G2 = Gprime.bottomRows(C.rows() - k);
    auto G2_LU = G2.fullPivLu();
    MatrixD G2_U = (G2_LU.matrixLU().template triangularView<Upper>());

    // Now we have a U matrix with zeros in the bottom left corner, while algorithm
    // wants them (at this point) in the upper left.  Reverse the rows to accomodate
    typedef PermutationMatrix<Dynamic, Dynamic, std::size_t> PermutationD;
    PermutationD reverse_rows;                // order of rows is completely reversed
    reverse_rows.setIdentity(G2.rows());      // start with null permutation
    for (std::size_t i = 0; i < (G2.rows() / 2); ++i) {
        reverse_rows.applyTranspositionOnTheRight(i, (G2.rows()-1) - i);
    }
    // permute columns of the full matrix according to LU result
    MatrixD Gnew = Gprime * G2_LU.permutationQ();
    // insert already-permuted rows that came from LU, but in reverse order
    Gnew.block(k, 0, Gprime.rows() - k, Gprime.cols()) =
        reverse_rows * G2_U;

    // Step 4: "Carry out the same row operations in the B matrix"
    // Note: not necessary to do it for C, because all coefficients are zero in those rows
    MatrixD G2_L = G2_LU.matrixLU().leftCols(G2.rows()).template triangularView<UnitLower>();
    MatrixD Bnew = Bprime;
    Bnew.block(k, 0, Bnew.rows() - k, Bnew.cols()) =
        reverse_rows * G2_L.inverse() * G2_LU.permutationP() * Bprime.bottomRows(Bprime.rows() - k);

    // Step 5: "Interchange the columns in the G, C, and D matrices... such that G22 is non-singular"
    // This comes down to applying both the column permutation from the LU decomposition of G2
    // and the column permutation that moves the resulting zeros into position in G22
    PermutationD exchange_columns;
    exchange_columns.setIdentity(G2.cols());
    for (std::size_t i = 0; i < (G2.cols() - k); ++i) {
        exchange_columns.applyTranspositionOnTheRight(i, (G2.cols()-1) - i);
    }
    Gnew = Gnew * exchange_columns;   // already has LU column pivots applied
    MatrixD Cnew = Cprime * G2_LU.permutationQ() * exchange_columns;
    Matrix<Float, ocount, Dynamic> Dnew = Dprime * G2_LU.permutationQ() * exchange_columns;

    // Step 6: compute reduced matrices using equations given in paper
    MatrixD G11 = Gnew.topLeftCorner(k, k);
    MatrixD G12 = Gnew.topRightCorner(k, Gnew.rows() - k);
    MatrixD G21 = Gnew.bottomLeftCorner(Gnew.rows() - k, k);
    MatrixD G22 = Gnew.bottomRightCorner(Gnew.rows() - k, Gnew.rows() - k);
    MatrixD C11 = Cnew.topLeftCorner(k, k);
    MatrixD C12 = Cnew.topRightCorner(k, Cnew.rows() - k);
    MatrixD B1  = Bnew.topRows(k);
    MatrixD B2  = Bnew.bottomRows(Bnew.rows() - k);
    MatrixD D1  = Dnew.leftCols(k);
    MatrixD D2  = Dnew.rightCols(Dnew.cols() - k);

    assert(!isSingular(G22));
    auto    G22_LU = G22.fullPivLu();

    MatrixD Gfinal  = G11 - G12 * G22_LU.solve(G21);
    MatrixD Cfinal  = C11 - C12 * G22_LU.solve(G21);
    MatrixD Bfinal  = B1  - G12 * G22_LU.solve(B2);
    MatrixD Bfinal2 =     - C12 * G22_LU.solve(B2);
    MatrixD Dfinal  = D1  - D2  * G22_LU.solve(G21);

    MatrixD Efinal  = D2 * G22_LU.solve(B2);

    // Checks
    assert(Efinal.isZero());         // assuming no feedthrough
    assert(!isSingular(Cfinal));     // assuming one pass is enough

    return std::make_tuple(Gfinal, Cfinal, Bfinal,
                           Dfinal.transpose());  // for PRIMA compatibility

}


}}}  // close namespaces

#endif  // EDASKEL_ANALYSIS_MNA_HPP

