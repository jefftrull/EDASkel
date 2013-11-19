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


// take a circuit's linear system description in G, C, B, L form and compress it so
// the resulting C array is non-singular.  Operation depends on runtime data, so
// output array dimensions are "Dynamic"
template<int icount, int ocount, int scount, typename Float = double>
std::tuple<Matrix<Float, Dynamic, Dynamic>,   // G result
           Matrix<Float, Dynamic, Dynamic>,   // C result
           Matrix<Float, Dynamic, icount>,    // B result
           Matrix<Float, Dynamic, ocount> >   // L result
regularize(Matrix<Float, scount, scount> const & G,
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
   typedef Matrix<Float, scount, scount> eqn_matrix_t;
   eqn_matrix_t permut = eqn_matrix_t::Identity();   // start with null permutation
   std::size_t i, j;
   for (i = 0, j=(scount-1); i < j;) {
      // loop invariant: rows > j are all zero; rows < i are not
      while ((i < scount) && !zero_rows(i)) ++i;
      while ((j > 0) && zero_rows(j)) --j;
      if (i < j) {
         // exchange rows i and j via the permutation vector
         permut(i, i) = 0; permut(j, j) = 0;
         permut(i, j) = 1; permut(j, i) = 1;
         ++i; --j;
      }
   }

   // 2. Apply permutation to MNA matrices
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

   auto G22QR = G22.fullPivLu();
   MatrixD G22invG21 = G22QR.solve(G21);
   auto G22invB2 = G22QR.solve(B2);
   MatrixD Gred = G11 - G12 * G22invG21;

   Matrix<Float, Dynamic, ocount> Lred = (L1.transpose() - L2.transpose() * G22invG21).transpose();
   Matrix<Float, Dynamic, icount> Bred = B1 - G12 * G22invB2;

   return std::make_tuple(Gred, Cred, Bred, Lred);
}

}}}  // close namespaces
         
#endif  // EDASKEL_ANALYSIS_MNA_HPP

