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
   // A singular matrix has at least one zero eigenvalue -
   // in theory, at least... but due to machine precision we can have "nearly singular"
   // matrices that misbehave.  Comparing rank instead is safer because it uses thresholds
   // for near-zero values.

   assert(m.rows() == m.cols());   // singularity has no meaning for a non-square matrix
   return (m.fullPivLu().rank() != m.rows());

}

// Calculate moments of given system in MNA form
template<int icount, int ocount, int scount, typename Float = double>
std::vector<Matrix<Float, ocount, icount>>
moments(Matrix<Float, scount, scount> const & G,
        Matrix<Float, scount, scount> const & C,
        Matrix<Float, scount, icount> const & B,
        Matrix<Float, scount, ocount> const & L,
        Matrix<Float, ocount, icount> const & E,
        size_t count) {
    std::vector<Matrix<Float, ocount, icount>> result;

    using namespace EDASkel::analysis::mna;
    assert(!isSingular(G));
    auto G_QR = G.fullPivHouseholderQr();
    Matrix<Float, scount, scount> A = -G_QR.solve(C);
    Matrix<Float, scount, icount> R = G_QR.solve(B);

    result.push_back(L.transpose() * R + E);   // incorporate feedthrough into first moment
    Matrix<Float, scount, scount> AtotheI = A;
    for (size_t i = 1; i < count; ++i) {
        result.push_back(L.transpose() * AtotheI * R);
        AtotheI = A * AtotheI;
    }

    return result;
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
   std::size_t nonzero_count = C.rows() - zero_count;

   // 1. Generate permutation matrix to move zero rows to the bottom
   PermutationMatrix<scount, scount, std::size_t> permut;
   permut.setIdentity(C.rows());      // start with null permutation
   std::size_t i, j;
   for (i = 0, j=(C.rows()-1); i < j;) {
      // loop invariant: rows > j are all zero; rows < i are not
      while ((i < C.rows()) && !zero_rows(i)) ++i;
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

   // This approach presumes no feedthrough (input-to-output) term
   MatrixD D = L2.transpose() * G22invB2;
   assert(D.isZero());

   return std::make_tuple(Gred, Cred, Bred, Lred);
}

template<int icount, int ocount, int scount, typename Float = double>
std::tuple<Matrix<Float, Dynamic, Dynamic>,   // G result
           Matrix<Float, Dynamic, Dynamic>,   // C result
           Matrix<Float, Dynamic, icount>,    // B result
           Matrix<Float, Dynamic, ocount>,    // D result
           Matrix<Float, ocount, icount> >    // E result (feedthrough)
regularize_natarajan(Matrix<Float, scount, scount> const & G,
           Matrix<Float, scount, scount> const & C,
           Matrix<Float, scount, icount> const & B,
           Matrix<Float, scount, ocount> const & D) {

    // Implements the algorithm in [Natarajan]
    // Circuits, Devices and Systems, IEE Proceedings G, June 1991

    typedef Matrix<Float, Dynamic, Dynamic> MatrixD;

    // Step 1: put C into "Row Echelon" form by performing LU factorization
    auto lu = C.fullPivLu();
    auto k = lu.rank();
    if (k == C.rows()) {
        // C is already non-singular
        Matrix<Float, ocount, icount>   E = Matrix<Float, ocount, icount>::Zero();
        return std::make_tuple(G, C, B, D, E);
    }

    MatrixD U = lu.matrixLU().template triangularView<Upper>();
    MatrixD L = lu.matrixLU().template triangularView<UnitLower>();

    // Step 2: "perform the same elementary operations on G and B"
    // given that C = P.inverse() * L * U * Q.inverse()
    // (from source it seems that permutationP/Q is inverse)
    // then to get the new G we reverse those operations:
    auto Cprime = U;   // note we may have small non-zero values in bottom rows, but they will be ignored
    auto P = lu.permutationP();
    auto Q = lu.permutationQ();

    assert(!isSingular(L));
    MatrixD Gprime = L.fullPivLu().solve(P * G * Q);                   // rows and columns
    Matrix<Float, scount, icount> Bprime = L.fullPivLu().solve(P * B); // rows only

    // The D input is like L in PRIMA but this algorithm uses the transpose
    Matrix<Float, ocount, scount> Dprime = D.transpose() * Q;          // columns only

    // Step 3: "Convert [G21 G22] matrix into row echelon form starting from the last row"
    MatrixD Cnew, Gnew, Bnew, Dnew;

    if (Cprime.rows() == (k+1)) {
        // if G22 is only a single row, there is no point attempting to decompose it
        Cnew = Cprime; Gnew = Gprime; Bnew = Bprime; Dnew = Dprime;
    } else {
        // decompose the bottom rows
        // Upon close review of the first example in the paper, the author is not only
        // converting from the last row, but also *from the last column*, i.e., he
        // performs a standard gaussian elimination on the matrix rotated 180 degrees

        // Plan of attack: reverse G2, perform LU decomposition, reverse result, reverse permutations
        MatrixD G2R = Gprime.bottomRows(C.rows() - k).reverse();

        auto G2R_LU = G2R.fullPivLu();
        MatrixD G2R_U = (G2R_LU.matrixLU().template triangularView<Upper>());

        // Since the order of the rows is irrelevant, I'll perform the decomposition, then
        // combine reversing the rows with the row reordering the LU decomposition produces

        MatrixD exchange_columns = G2R_LU.permutationQ();
        Gnew = Gprime * exchange_columns.reverse();

        // insert already-permuted rows that came from LU, but in reverse order
        Gnew.block(k, 0, Gprime.rows() - k, Gprime.cols()) = G2R_U.reverse();

        // Step 4: "Carry out the same row operations in the B matrix"
        // Note: not necessary to do it for C, because all coefficients are zero in those rows
        // 4.1 reverse the rows in B2
        typedef PermutationMatrix<Dynamic, Dynamic, std::size_t> PermutationD;
        PermutationD reverse_rows;                // order of rows is completely reversed
        reverse_rows.setIdentity(G2R.rows());     // start with null permutation
        for (std::size_t i = 0; i < (G2R.rows() / 2); ++i) {
            reverse_rows.applyTranspositionOnTheRight(i, (G2R.rows()-1) - i);
        }
        MatrixD B2R = reverse_rows * Bprime.bottomRows(Bprime.rows() - k);

        // 4.2 extract and apply L operation from reversed G2
        MatrixD G2R_L = G2R_LU.matrixLU().leftCols(G2R.rows()).template triangularView<UnitLower>();
        Bnew = Bprime;
        Bnew.block(k, 0, Bnew.rows() - k, Bnew.cols()) =
            reverse_rows.transpose() * G2R_L.fullPivLu().solve(G2R_LU.permutationP() * B2R);

        // Step 5: "Interchange the columns in the G, C, and D matrices... such that G22 is non-singular"
        // Since we have done a full pivot factorization of G2 I assume G22 is already non-singular,
        // so the only thing left to do is reorder the C and D matrices according to the G2 factorization
        Cnew = Cprime * exchange_columns.reverse();
        Dnew = Dprime * exchange_columns.reverse();
    }

    // Step 6: compute reduced matrices using equations given in paper
    MatrixD G11 = Gnew.topLeftCorner(k, k);
    MatrixD G12 = Gnew.topRightCorner(k, Gnew.rows() - k);
    MatrixD G21 = Gnew.bottomLeftCorner(Gnew.rows() - k, k);
    MatrixD G22 = Gnew.bottomRightCorner(Gnew.rows() - k, Gnew.rows() - k);
    MatrixD C11 = Cnew.topLeftCorner(k, k);
    MatrixD C12 = Cnew.topRightCorner(k, Cnew.rows() - k);
    MatrixD B01 = Bnew.topRows(k);
    MatrixD B02 = Bnew.bottomRows(Bnew.rows() - k);
    MatrixD D01 = Dnew.leftCols(k);
    MatrixD D02 = Dnew.rightCols(Dnew.cols() - k);

    assert(!isSingular(G22));
    auto    G22_LU = G22.fullPivLu();

    MatrixD Gfinal  = G11 - G12 * G22_LU.solve(G21);
    MatrixD Cfinal  = C11 - C12 * G22_LU.solve(G21);
    Matrix<double, Dynamic, icount> B1
                    = B01 - G12 * G22_LU.solve(B02);
    MatrixD B2      =     - C12 * G22_LU.solve(B02);
    Matrix<double, ocount, Dynamic> Dfinal
                   = D01 - D02 * G22_LU.solve(G21);

    Matrix<double, ocount, icount> E1
                    =       D02 * G22_LU.solve(B02);

    assert(!isSingular(Cfinal));   // not iterating yet b/c we need to combine results

    // Now apply a transformation suggested by Chen (TCAD July 2012) to eliminate
    // the input derivative term (B2)
    auto CfinalQR = Cfinal.fullPivHouseholderQr();
    Matrix<double, Dynamic, icount> Bfinal      = B1 - Gfinal * CfinalQR.solve(B2);
    // This change of variable creates an additional feedthrough term
    Matrix<double, ocount, icount>  feedthrough =      Dfinal * CfinalQR.solve(B2);

    return std::make_tuple(Gfinal, Cfinal, Bfinal,
                           Dfinal.transpose(),  // for PRIMA compatibility
                           E1 + feedthrough);
}


}}}  // close namespaces

#endif  // EDASKEL_ANALYSIS_MNA_HPP

