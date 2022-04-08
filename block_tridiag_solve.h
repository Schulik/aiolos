
#ifndef _BLOCK_TRIDIAG_SOLVE_H
#define _BLOCK_TRIDIAG_SOLVE_H

#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>


class BlockTriDiagSolverBase
{
public:
  virtual void factor_matrix(double *l, double *d, double *u) = 0;
  virtual void solve(double *rhs, double *x) = 0;

  virtual ~BlockTriDiagSolverBase(){};
};

/* Solver a tri-diagonal systems of equations where each
 * of the blocks have the same size (N x N).
 * 
 *
 */
template <int fixed_block_size = Eigen::Dynamic>
class BlockTriDiagSolver
    : public BlockTriDiagSolverBase
{
  // Inernal vector/matrix views
  typedef Eigen::Map<Eigen::Matrix<double, fixed_block_size, 1>> vec_t;
  typedef Eigen::Map<Eigen::Matrix<double, 
                      fixed_block_size, fixed_block_size,
                      Eigen::RowMajor>> mat_t;
public:
  BlockTriDiagSolver() {} ;
  BlockTriDiagSolver(int num_blocks)
      : _sLU(num_blocks),
        _l(num_blocks * fixed_block_size * fixed_block_size),
        _u(num_blocks * fixed_block_size * fixed_block_size),
        _n_blocks(num_blocks), _block_size(fixed_block_size){};

  BlockTriDiagSolver(int num_blocks, int block_size)
      : _sLU(num_blocks),
        _l(num_blocks * block_size * block_size),
        _u(num_blocks * block_size * block_size),
        _n_blocks(num_blocks), _block_size(block_size)
  {
    if (fixed_block_size != Eigen::Dynamic)
      throw std::invalid_argument("block_size parameter should not be passed "
                                  "to the constructor unless dynamic memory is "
                                  "used.");
  };

  /* Pre-compute the matrix factorization for later solution */
  void factor_matrix(double *l, double *d, double *u)
  {
    int step = _block_size * _block_size;

    Eigen::Matrix<double, fixed_block_size, fixed_block_size> S_i =
        mat_t(d, _block_size, _block_size);

    // Compute the LU decomposition
    _sLU[0].compute(S_i);

    for (int i = 1; i < _n_blocks; i++)
    {
      // Load the required matrices
      mat_t l_0(l + i * step, _block_size, _block_size);
      mat_t d_0(d + i * step, _block_size, _block_size);
      mat_t u_m(u + (i - 1) * step, _block_size, _block_size);

      // Update d_0 to give the next denominator:
      S_i.noalias() = d_0 - l_0 * _sLU[i-1].solve(u_m);
      // std::cout << i << ":\n" << d_0 << "\n\n" << S_i << std::endl ;

    // Compute the LU decomposition
      _sLU[i].compute(S_i);
    }

    // Save the l/u matrices for later
    for (int i = 0; i < _n_blocks * step; i++)
    {
      _l[i] = l[i];
      _u[i] = u[i];
    }
  }

  void solve(double *rhs, double *x)
  {
    int step = _block_size * _block_size;

    double *l = &_l.front();
    double *u = &_u.front();

    // Do the forward-pass:
    //   Replace the vectors on the RHS with the vectors
    //      A_i = [S_i]^-1 * (rhs_i - l_i * A_i-1),
    //   Note l_0 = 0.
    vec_t rhs_0(rhs, _block_size, 1) ;
    rhs_0 = _sLU[0].solve(rhs_0).eval() ;

    for (int i = 1; i < _n_blocks; i++)
    {
      // Load the required matrix/vectors
      vec_t rhs_m(rhs + (i-1)*_block_size, _block_size, 1);
      vec_t rhs_0(rhs +   i  *_block_size, _block_size, 1);
      mat_t l_0(l + i*step, _block_size, _block_size);

      rhs_0 = _sLU[i].solve(rhs_0 - l_0*rhs_m).eval() ;
    }

    // Now do the back-subsitution for x@
    //   x_i = A_i - [S_i]^-1 * u_i * x_i+1
    vec_t x_n(x + (_n_blocks - 1) * _block_size, _block_size, 1);
    vec_t rhs_n(rhs + (_n_blocks - 1) * _block_size, _block_size, 1);

    x_n = rhs_n ;

    for (int i = _n_blocks - 2; i >= 0; i--)
    {
      vec_t x_p(x + (i + 1) * _block_size, _block_size, 1);
      vec_t x_0(x + i * _block_size, _block_size, 1);
      vec_t rhs_0(rhs + i * _block_size, _block_size, 1);

      mat_t u_0(u + i * step, _block_size, _block_size);

      x_0.noalias() = rhs_0 - _sLU[i].solve(u_0 * x_p);
    }
  }

private:
  typedef Eigen::PartialPivLU<Eigen::Matrix<double, fixed_block_size,
                                            fixed_block_size>>
      LU;
  std::vector<LU> _sLU;
  std::vector<double> _l, _u;

  int _n_blocks, _block_size;
};

/* Solver a tri-diagonal systems of equations where each
 * of the blocks have the same size (N x N), using Eigens
 * Sparse matrix solvers
 *
 */
template <int fixed_block_size = Eigen::Dynamic>
class BlockTriDiagSolverSparse
    : public BlockTriDiagSolverBase
{
  // Inernal vector/matrix views
  typedef Eigen::Map<Eigen::Matrix<double, fixed_block_size, 1>> vec_t;

public:
  BlockTriDiagSolverSparse() {} ;
  BlockTriDiagSolverSparse(int num_blocks)
      : _mat(num_blocks*fixed_block_size, num_blocks*fixed_block_size),
        _n_blocks(num_blocks), _block_size(fixed_block_size) 
  {
    _mat.reserve(Eigen::VectorXi::Constant(num_blocks*fixed_block_size, 3*fixed_block_size)) ;
  } ;

  BlockTriDiagSolverSparse(int num_blocks, int block_size)
      : _mat(num_blocks*block_size, num_blocks*block_size),
        _n_blocks(num_blocks), _block_size(block_size)
  {
    _mat.reserve(Eigen::VectorXi::Constant(num_blocks*block_size, 3*block_size)) ;

    if (fixed_block_size != Eigen::Dynamic)
      throw std::invalid_argument("block_size parameter should not be passed "
                                  "to the constructor unless dynamic memory is "
                                  "used.");
  };

  /* Pre-compute the matrix factorization for later solution */
  void factor_matrix(double *l, double *d, double *u)
  {
    // Setup the matrix
    _mat.setZero() ;
    for (int i=0; i <_n_blocks; i++) {
      for (int j=0; j < _block_size; j++) {
        for(int k=0; k < _block_size; k++) {
          int idx = i*_block_size*_block_size + j*_block_size + k ;
            
          if (i > 0)
            _mat.insert(i*_block_size+j, (i-1)*_block_size+k) = l[idx] ;
          _mat.insert  (i*_block_size+j,   i  *_block_size+k) = d[idx] ;
          if (i+1 < _n_blocks)
            _mat.insert(i*_block_size+j, (i+1)*_block_size+k) = u[idx] ;
        } 
      }  
    }

    _solver.compute(_mat) ;

    if (_solver.info() != Eigen::Success) 
        throw std::runtime_error("Factorization failed") ;
  }

  void solve(double* _rhs, double* _x)
  {
    vec_t rhs(_rhs, _block_size*_n_blocks, 1) ;
    vec_t x(_x, _block_size*_n_blocks, 1) ;

    x = _solver.solve(rhs) ;
  } ;


  BlockTriDiagSolverSparse<fixed_block_size>& operator=(BlockTriDiagSolverSparse<fixed_block_size>&& o) {
    _mat = std::move(o._mat) ;
    _n_blocks = o._n_blocks ;
    _block_size = o._block_size ;

    return *this ;
  }

private:
  Eigen::SparseMatrix<double> _mat ;
  
  Eigen::SparseLU<decltype(_mat), Eigen::COLAMDOrdering<int>> _solver ;
  

  int _n_blocks, _block_size;
};


#endif//_BLOCK_TRIDIAG_SOLVE_H