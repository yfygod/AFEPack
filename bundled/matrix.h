/**
 * @file matrix.h
 * @brief 矩阵类，继承自Eigen中的Matrix
 * @author Yang Fanyi <yfy__92@126.com>
 * @version 
 */

#ifndef __MATRIX_H_
#define __MATRIX_H_

#include "vector.h"

template <typename T, int R, int C, int Major = RowMajor>
class MatrixType : public Eigen::Matrix<T, R, C, Major> {
  typedef Eigen::Matrix<T, R, C, Major> base_t;
  typedef typename base_t::Index index_t; 

public:
  enum {
    RowsAtCompileTime = base_t::RowsAtCompileTime,
    ColsAtCompileTime = base_t::ColsAtCompileTime,
    SIZE = (R > C) ? R : C,
    SIZE0 = R, SIZE1 = C, 
  };

  // constructors:

public:
  MatrixType() : base_t() {}

#ifdef EIGEN_HAVE_RVALUE_REFERENCES
  MatrixType(MatrixType&& other) : base_t(std::move(other)){}
  MatrixType& operator = (MatrixType&& other) {
    other.swap(*this); return *this;}
#endif

  explicit MatrixType(const index_t l) : base_t(l) {}

#ifndef EIGEN_PARSED_BY_DOXYGEN
  template <typename T0, typename T1>
  inline MatrixType(const T0& x, const T1& y) : base_t(x, y){}
#else
  MatrixType(index_t r, index_t c) : base_t(r, c) {}
  MatrixType(const T& x, const T& y) : base_t(x, y){}
#endif

  MatrixType(const T& x, const T& y, const T& z) : 
    base_t(x, y, z) {}
  MatrixType(const T& x, const T& y, const T& z, const T& w) : 
    base_t(x, y, z, w) {}

  template <typename OtherDerived>
  inline MatrixType(const Eigen::MatrixBase<OtherDerived>& other) :
    base_t(other) {}

  inline MatrixType(const MatrixType& other) : base_t(other) {}

  template <typename OtherDerived>
  inline MatrixType(const Eigen::ReturnByValue<OtherDerived>& other) : 
    base_t(other) {}

  template <typename OtherDerived>
  inline MatrixType(const Eigen::EigenBase<OtherDerived>& other) : 
    base_t(other) {}

  inline MatrixType& operator = (const MatrixType& other) {
    base_t::operator = (other);
    return *this;
  }

  template <typename OtherDerived>
  inline MatrixType& operator = (const Eigen::MatrixBase<OtherDerived>& other) {
    base_t::operator = (other);
    return *this;
  }

  template <typename OtherDerived>
  inline MatrixType& operator = (const Eigen::EigenBase<OtherDerived>& other) {
    base_t::operator = (other);
    return *this;
  }

  template <typename OtherDerived>
  inline MatrixType& operator = (const Eigen::ReturnByValue<OtherDerived>& other){
    base_t::operator = (other);
    return *this;
  }

  template <typename OtherDerived>
  inline explicit MatrixType(const Eigen::RotationBase<OtherDerived, 
      ColsAtCompileTime>& r) : base_t(r){}

  template <typename OtherDerived>
  inline MatrixType& operator = (const Eigen::RotationBase<OtherDerived,
      ColsAtCompileTime>& r) {
    base_t::operator = (r);
    return *this;
  }

  inline void set_zero() { base_t::setZero(); }

  inline void set_constant(const T t) { base_t::setConstant(t); }
};

typedef MatrixType<double, 1, 1> Matrix1d; typedef MatrixType<double, 2, 2> Matrix2d;
typedef MatrixType<double, 3, 3> Matrix3d; typedef MatrixType<double, 4, 4> Matrix4d;
typedef MatrixType<double, 5, 5> Matrix5d; typedef MatrixType<double, 6, 6> Matrix6d;
typedef MatrixType<double, 7, 7> Matrix7d; typedef MatrixType<double, 8, 8> Matrix8d;

typedef MatrixType<double, Dynamic, Dynamic> Matrix;

// This class is required by AFEPack.
// This class substitues the class dealII::FullMatrix

template <typename T> 
class FullMatrix : public MatrixType<T, Dynamic, Dynamic> {
  typedef MatrixType<T, Dynamic, Dynamic> base_t;
  typedef typename Eigen::Matrix<T, Dynamic, Dynamic>::Index index_t;

public:
  using base_t::base_t; 

  // the matrix operations below are the same as dealII::FullMatrix

  template <typename V0, typename V1> 
  void vmult(V0& v0, const V1& v1) const {
    AssertExp2(base_t::rows() == v0.size(), base_t::rows(), v0.size());
    AssertExp2(base_t::cols() == v1.size(), base_t::cols(), v1.size());
    v0 = (*this)*v1;
  }

  void reinit(const index_t m, const index_t n) {
    base_t::resize(m, n); 
    base_t::setZero();
  }

  void gauss_jordan() {
    AssertExp2(base_t::rows() == base_t::cols(), base_t::rows(), base_t::cols());
    Eigen::FullPivLU<Eigen::MatrixXd> full_lu(*this);

    FullMatrix<T> tmp(base_t::rows(), base_t::cols());
    tmp = full_lu.inverse(); (*this) = tmp;
  }

};


#endif
