/**
 * @file vector.h
 * @brief 向量类，继承自Eigen中的Matrix
 * @author Yang Fanyi <yfy__92@126.com>
 * @version 
 */

#ifndef __VECTOR_H_
#define __VECTOR_H_

#ifndef DEBUG
#define EIGEN_NO_DEBUG
#endif

#include "exception.h"
#include "Eigen/Dense"

#include <fstream>

enum {
  Dynamic = Eigen::Dynamic, 
  ColMajor = Eigen::ColMajor, 
  RowMajor = Eigen::RowMajor,
};

template <typename T, int S> 
class VectorType : public Eigen::Matrix<T, S, 1> {
  public:
  typedef Eigen::Matrix<T, S, 1> base_t;
  typedef typename base_t::Index Index;

public:
  enum {isDynamic = (S > Dynamic) ? false : true, 
    SIZE = (S > Dynamic) ? S : Dynamic,
    ColsAtCompileTime = base_t::ColsAtCompileTime,
  };

  // constructors
public:
  VectorType() : base_t() {}

  explicit VectorType(const Index l) : base_t(l) {
    if(isDynamic) base_t::setZero();
  }

#ifndef EIGEN_PARSED_BY_DOXYGEN
  template <typename T0, typename T1>
  VectorType(const T0& x, const T1& y) : base_t(x, y){}
#else
  VectorType(Index rows, Index cols) : base_t(rows, cols) {}
  VectorType(const T& x, const T& y) : base_t(x, y) {}
#endif

  VectorType(const T& x, const T& y, const T& z) : 
    base_t(x, y, z) {}

  VectorType(const T& x, const T& y, const T& z, const T& w) : 
    base_t(x, y, z, w) {}

  template <typename OtherDerived> inline
  VectorType(const Eigen::MatrixBase<OtherDerived>& other) : 
    base_t(other) {}

  inline
  VectorType(const VectorType& other) : base_t(other) {}

  template <typename OtherDerived> inline
  VectorType(const Eigen::EigenBase<OtherDerived>& other) : 
    base_t(other) {}

  template <typename OtherDerived> inline
  explicit VectorType(const Eigen::RotationBase<OtherDerived, 
      ColsAtCompileTime>& r) : base_t(r){}

  //default operator

  inline VectorType& operator = (const VectorType& other) {
    base_t::operator = (other); 
    return *this;
  }

  template <typename OtherDerived> inline
  VectorType& operator = (const Eigen::MatrixBase<OtherDerived>& other) {
    base_t::operator = (other);
    return *this;
  }

  inline void set_zero() { base_t::setZero(); }

  inline void set_constant(const T t) { base_t::setConstant(t); }

  void pritn(std::ostream& os) const {
    os << (*this) << "\n";
  }

  void print(const std::string& file, int precision = 13) const {
    std::ofstream fout(file); fout.precision(precision);
    fout << (*this); fout.close();
  }
};

typedef VectorType<double, 1> Vector1d; typedef VectorType<double, 2> Vector2d;
typedef VectorType<double, 3> Vector3d; typedef VectorType<double, 4> Vector4d;
typedef VectorType<double, 5> Vector5d; typedef VectorType<double, 6> Vector6d;
typedef VectorType<double, 7> Vector7d; typedef VectorType<double, 8> Vector8d;

// This class is required by AFEPack.
// This class substitues the class dealII::Vector

template <typename T> 
class Vector : public VectorType<T, Dynamic> {
  typedef VectorType<T, Dynamic> base_t;
  typedef typename Eigen::Matrix<T, Dynamic, 1>::Index index_t;
public:
  using base_t::base_t; 

  typedef T* iterator;
  typedef const T* const_iterator; 
  typedef T& reference; 
  typedef const T& const_reference; 

  Vector& operator = (const Vector& other) {
    base_t::operator = (other);
    return *this;
  }

  void operator = (const T t) {
    base_t::set_constant(t); 
  }

  void reinit(const index_t l) {
    base_t::resize(l);
    base_t::set_zero(); 
  }

  template <typename stream>
  void block_read(stream& is) { 
    for(int i = 0; i < base_t::size(); ++i)
      is >> base_t::operator[](i);
  }

  template <typename stream>
  void block_write(stream& os) const { 
    for(int i = 0; i < base_t::size(); ++i)
      os << base_t::operator[](i);
  }

  // the vector operations below are the same as dealII::Vector

  void add(const double r, const Vector& other) {
    AssertExp2(base_t::size() == other.size(), base_t::size(), other.size());
    for(int i = 0; i < base_t::size(); ++i)
      (*this)(i) += (r*other(i));
  }

  void addv(const double a, const Vector& other, const double b) {
    AssertExp2(base_t::size() == other.size(), base_t::size(), other.size());
    T tmp; 
    for(int i = 0; i < base_t::size(); ++i) {
      tmp = (*this)(i);
      (*this)(i) = (a*tmp + b*other(i));
    }
  }

  void sadd(const double r, const Vector& other) {
    (*this) *= r; 
    (*this) += other;
  }

  void sadd(const double a, const double b, const Vector& other) {
    (*this) *= a; (*this) += (b*other);
  }

  void sadd(const double s, 
      const double a,
      const Vector& V,
      const double b,
      const Vector& W) {
    (*this) *= s; (*this) += (a*V) + (b*W);
  }

  double l1_norm() const { return (*this).template lpNorm<1>();}
  double l2_norm() const { return (*this).template lpNorm<2>();}

  bool all_zero() const {
    return base_t::isZero(1.e-15);
  }

  iterator begin() { return base_t::data(); }
  const_iterator begin() const { return base_t::data(); }

  iterator end() { return base_t::data() + base_t::size();}
  const_iterator end() const { return base_t::data() + base_t::size();}

};


#endif
