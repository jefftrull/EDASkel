// An implementation of the PRIMA model order reduction algorithm using Eigen
// part of EDASkel, a sample EDA app
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
#if !defined(EDASKEL_ANALYSIS_PRIMA_HPP)
#define EDASKEL_ANALYSIS_PRIMA_HPP

#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>

namespace EDASkel {
namespace analysis {

// PRIMA
// This is a famous model reduction technique invented around 1997/8 at CMU
// I am generally following the treatment in Odabasioglu, IEEE TCAD, August 1998
// They base their work on the prior Block Arnoldi algorithm, a helpful
// explanation of which you can find in their 6th citation:
// D. L. Boley, “Krylov space methods on state-space control models”
      
template<typename Float>
Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>
Prima(Eigen::SparseMatrix<Float> const& C,   // derivative conductance terms
      Eigen::SparseMatrix<Float> const& G,   // conductance
      Eigen::SparseMatrix<Float> const& B,   // input
      Eigen::SparseMatrix<Float> const& L,   // output
      std::size_t q) {                       // desired state variables
  // assert preconditions
  assert(C.rows() == C.cols());     // input matrices are square
  assert(G.rows() == G.cols());
  assert(C.rows() == G.rows());     // input matrices are of the same size
  assert(B.rows() == G.rows());
  assert(L.rows() == G.rows());
  size_t N = B.cols();
  assert(N < C.rows());             // must have more state variables than ports
  assert(q < C.rows());             // desired state count must be less than current number

  // unchecked precondition: the state variables associated with the ports must be the last N

  using namespace Eigen;
  using namespace std;

  // Step 1 of PRIMA creates the B and L matrices, and is performed by the caller.

  // Step 2: Solve GR = B for R
  SparseLU<SparseMatrix<Float>, COLAMDOrdering<int> > G_QR(G);
  assert(G_QR.info() == Success);
  SparseMatrix<Float> R = G_QR.solve(B);
  assert(G_QR.info() == Success);

  // Step 3: Set X[0] to the orthonormal basis of R as determined by QR factorization
  typedef Matrix<Float, Dynamic, Dynamic> MatrixXX;
  // The various X matrices are stored in a std::vector.  Eigen requires us to use a special
  // allocator to retain alignment:
  typedef aligned_allocator<MatrixXX> AllocatorXX;
  typedef vector<MatrixXX, AllocatorXX> MatrixXXList;
  
  SparseQR<SparseMatrix<Float>, COLAMDOrdering<int> > R_QR(R);
  assert(R_QR.info() == Success);
  SparseMatrix<Float> rQ;
  rQ = R_QR.matrixQ();
  MatrixXXList X(1, rQ.leftCols(R_QR.rank()));

  // Step 4: Set n = floor(q/N)+1 if q/N is not an integer, and q/N otherwise
  size_t n = (q % N) ? (q/N + 1) : (q/N);

  // Step 5: Block Arnoldi (see Boley for detailed explanation)
  // In some texts this is called "band Arnoldi".
  // Boley and PRIMA paper use X with both subscripts and superscripts
  // to indicate the outer (subscript) and inner (superscript) loops
  // I have used X[] for the outer, Xk[] for the inner

  // pre-calculate G^-1*C for use in inner loop
  SparseMatrix<Float> GinvC = G_QR.solve(C);

  for (size_t k = 1; k < n; ++k)
  {
    // because X[] will vary in number of columns, so will Xk[]
    MatrixXXList Xk(k+1);

    // Prima paper says:
    // set V = C * X[k-1]
    // solve G*X[k][0] = V for X[k][0]

    // X[k][0] = G^-1*V = G^-1*C*X[k-1] = (G^-1*C)*X[k-1]
    Xk[0] = GinvC*X[k-1];     // So Xk[0] = G^-1*C*X[k-1], i.e. A*X[k-1]
                              // Boley: "expand Krylov space"

    for (size_t j = 1; j <= k; ++j)   // "Modified Gram-Schmidt"
    {
      auto H = X[k-j].transpose() * Xk[j-1];   // H[k-j][k-1] per Boley

      // X[k][j] = X[k][j-1] - X[k-j]*H
      Xk[j] = Xk[j-1] - X[k-j] * H;
    }

    // set X[k] to the orthonormal basis of X[k][k] via QR factorization
    // per Boley the "R" produced is H[k][k-1]
    if (Xk[k].cols() == 1)
    {
      // a single column is automatically orthogonalized; just normalize
      X.push_back(Xk[k].normalized());
    } else {
      auto xkkQR = Xk[k].fullPivHouseholderQr();
      MatrixXX xkkQ = xkkQR.matrixQ();
      X.push_back(xkkQ.leftCols(xkkQR.rank()));
    }
  }

  // Step 6: Set Xfinal to the concatenation of X[0] to X[n-1],
  //         truncated to q columns
  size_t cols = accumulate(X.begin(), X.end(), 0,
                           [](size_t sum, MatrixXX const& m) { return sum + m.cols(); });
  cols = std::min(q, cols);  // truncate to q

  MatrixXX Xfinal(C.rows(), cols);
  size_t col = 0;
  for (size_t k = 0; (k <= n) && (col < cols); ++k)
  {
    // copy columns from X[k] to Xfinal
    for (int j = 0; (j < X[k].cols()) && (col < cols); ++j)
    {
      Xfinal.col(col++) = X[k].col(j);
    }
  }

  return Xfinal;

}
      
}}

#endif  // EDASKEL_ANALYSIS_PRIMA_HPP

