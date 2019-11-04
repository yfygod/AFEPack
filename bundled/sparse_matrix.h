/**
 * @file sparse_matrix.h
 * @brief SparseMatrix, it borrows from dealII::SparseMatrix
 * @author Yang Fanyi <yfy__92@126.com>
 * @version 
 * @date 2019-11-03
 */

#include <cstring>

#include "sparsity_pattern.h"

template <typename Scalar, typename Index = unsigned int>
class SparseMatrix {
public:
  typedef Index index_t;
  typedef Scalar scalar_t; 

  typedef scalar_t value_type;

  SparseMatrix() : cols(NULL), val(NULL), max_len(0) {}

  SparseMatrix(const SparseMatrix& other) : cols(NULL), val(NULL), 
                                            max_len(0) {} 

  SparseMatrix(const SparsityPatternBase<Index>& sp) {
    cols = NULL;
    val = NULL;
    max_len = 0;
    reinit(sp);
  }

  SparseMatrix& operator = (const SparseMatrix&) {
    return *this;
  }

  ~SparseMatrix() {
    cols = NULL; 
    if(val != NULL)
      delete val;
  }

public: 
  SparseMatrix& operator = (const Scalar);

  SparseMatrix& operator += (const SparseMatrix&);

  void reinit(const SparsityPatternBase<Index>&); 

  void clear(); 

  index_t n_actually_nonzero_elements(const double) const;

  void set(const index_t, const index_t, const scalar_t);

  void set_row(const index_t, const scalar_t);

  void add(const index_t, const index_t, const scalar_t);

  void set_zero(); 

  SparseMatrix& operator *= (const scalar_t);
  SparseMatrix& operator /= (const scalar_t);

  template <typename OUTVEC, typename INVEC>
  void vmult(OUTVEC&, const INVEC&) const;

  template <typename OUTVEC, typename INVEC>
  void Tvmult(OUTVEC&, const INVEC&) const;

  template <typename INVEC1, typename INVEC2>
  scalar_t vTmultv(const INVEC1&, const INVEC2&) const;

  template <typename VEC>
  double residual(const VEC&, const VEC&) const;

  void print(const std::string&, int = 13) const;

public:
  const SparsityPatternBase<Index>& get_sparsity_pattern() const {
    AssertCond(cols != NULL);
    return *cols;
  }

  bool empty() const {
    if(cols == NULL) return true;
    return cols->empty();
  }

  index_t m() const {
    if(cols == NULL) return 0;
    return cols->rows;
  }

  index_t n() const {
    if(cols == NULL) return 0;
    return cols->cols;
  }

  index_t n_nonzero_elements() const {
    if(cols == NULL) return 0;
    return cols->n_nonzero_elements();
  }

  index_t row_length(const index_t i) const {
    if(cols == NULL) return 0;
    return cols->row_length(i);
  }

  scalar_t operator() (const index_t i, const index_t j) const {
    Assert(cols != NULL, "should be initialized first");
    index_t index = cols->operator()(i, j);

    if(index != cols->invalid_entry)
      return val[index];
    else
      return 0;
  }

  scalar_t el(const index_t i, const index_t j) const {
    return this->operator()(i, j);
  }

  index_t max_entries_per_row() const {
    if(cols == NULL) return 0;
    return cols->max_entries_per_row();
  }

  const std::size_t * get_rowstart_indices() const {
    AssertCond(cols != NULL);
    return cols->get_rowstart_indices();
  }

  const index_t * get_column_numbers() const {
    AssertCond(cols != NULL);
    return cols->get_column_numbers(); 
  }

  scalar_t * get_nonzero_value() {
    return val;
  }

  const scalar_t * get_nonzero_value() const {
    return val;
  }

  scalar_t global_entry(const index_t i) const {
    Assert(cols != NULL, "not initialized");
    Assert(i < cols->n_nonzero_elements(), "exceed range");
    return val[i];
  }

  scalar_t& global_entry(const index_t i) {
    Assert(cols != NULL, "not initialized");
    Assert(i < cols->n_nonzero_elements(), "exceed range");
    return val[i];
  }

  scalar_t diag_element(const index_t i) const {
    Assert(cols != NULL, "not initialized");
    Assert(m() == n(), "not quadratic");
    Assert(i < m(), "exceed range");

    return val[cols->rowstart[i]];
  }

  scalar_t& diag_element(const index_t i) {
    Assert(cols != NULL, "not initialized");
    Assert(m() == n(), "not quadratic");
    Assert(i < m(), "exceed range");

    return val[cols->rowstart[i]];
  }


private:
  const SparsityPatternBase<Index> * cols; 
  scalar_t * val; 
  std::size_t max_len; 
};

template <typename ostream, typename Scalar, typename Index>
ostream& operator << (ostream& os, const SparseMatrix<Scalar, Index>& matrix) {
  const Scalar * val = matrix.get_nonzero_value(); 
  const std::size_t * rowstart = matrix.get_rowstart_indices();
  const Index * colnums = matrix.get_column_numbers(); 

  AssertCond(rowstart != NULL && colnums != NULL && val != NULL);

  for(Index i = 0; i < matrix.m(); ++i)
    for(std::size_t j = rowstart[i]; j < rowstart[i + 1]; ++j)
      os << '(' << i << ',' << colnums[j] << ')' << " " << val[j] << "\n";
    
  return (os << std::flush);
}

#include "sparse_matrix.templates.h"
