/**
 * @file sparsitypattern.h
 * @brief Sparse Pattern, it borrows from dealII::SparsityPattern
 * @author Yang Fanyi <yfy__92@126.com>
 * @version 
 * @date 2019-11-03
 */

#ifndef __SPARSITY_PATTERN_H_
#define __SPARSITY_PATTERN_H_

#include <vector>
#include <fstream>
#include <algorithm>

#include "exception.h"

template <typename Index = unsigned int>  
class SparsityPatternBase {
public:
  typedef Index index_t; 

  const index_t invalid_entry = static_cast<index_t>(-1);

  // constructors: 
public:
  SparsityPatternBase() : max_dim(0), 
                          max_vec_len(0), 
                          rowstart(NULL), 
                          colnums(NULL), 
                          compressed(false),
                          store_diagonal_first_in_row(false) {
    reinit(0, 0, 0);
  }

  // these two functions do not generate a copy form a given SparsityPatternBase

  SparsityPatternBase(const SparsityPatternBase&) : max_dim(0), 
                          max_vec_len(0), 
                          rowstart(NULL), 
                          colnums(NULL), 
                          compressed(false),
                          store_diagonal_first_in_row(false) {
    reinit(0, 0, 0);
  }

  SparsityPatternBase& operator = (const SparsityPatternBase&) {
    return *this;
  }

  ~SparsityPatternBase() {
    if(colnums != NULL) delete colnums;
    if(rowstart != NULL) delete rowstart; 
  } 

  SparsityPatternBase(const index_t m, const index_t max_per_row) : 
    max_dim(0),
    max_vec_len(0), 
    rowstart(NULL), 
    colnums(NULL), 
    compressed(false), 
    store_diagonal_first_in_row(true) {
      reinit(m, m, max_per_row);
  }

  SparsityPatternBase(const index_t m, const index_t n, 
      const index_t max_per_row) : 
    max_dim(0),
    max_vec_len(0), 
    rowstart(NULL), 
    colnums(NULL), 
    compressed(false), 
    store_diagonal_first_in_row(m == n) {
      reinit(m, n, max_per_row);
  }

  SparsityPatternBase(const index_t m, 
    const std::vector<index_t>& row_length) : 
    max_dim(0),
    max_vec_len(0), 
    rowstart(NULL), 
    colnums(NULL), 
    compressed(false), 
    store_diagonal_first_in_row(true) {
      reinit(m, m, row_length);
  }

  SparsityPatternBase(const index_t m, 
    const index_t n, 
    const std::vector<index_t>& row_length) : 
    max_dim(0),
    max_vec_len(0), 
    rowstart(NULL), 
    colnums(NULL), 
    compressed(false), 
    store_diagonal_first_in_row(m == n) {
      reinit(m, n, row_length);
  }

  SparsityPatternBase(const SparsityPatternBase&, 
      const index_t, const index_t);

public:
  void reinit(const index_t, const index_t, const index_t); 

  void reinit(const index_t, const index_t, const std::vector<index_t>&); 

  void add(const index_t, const index_t);

  void compress(); 

  index_t operator() (const index_t, const index_t) const; 

public:
  index_t n_rows() const { return rows;}
  index_t n_cols() const { return cols;}

  bool is_compressed() const { return compressed;}

  bool empty() const {
    if((rowstart == NULL) || (rows == 0) || (cols == 0)){
      return true;
    }
    return false;
  }

  index_t n_nonzero_elements() const { 
    AssertExc( (rowstart != NULL) && ( colnums != NULL),
        "the sparsity pattern is empty");
    AssertExc(compressed, "the sparsity pattern is uncompressed");
    return (rowstart[rows] - rowstart[0]);
  }

  index_t max_entries_per_row() const {
    if(!compressed)
      return max_row_length;

    index_t m = 0;
    for(index_t i = 1; i < rows + 1; ++i)
      m = std::max(m, static_cast<index_t>(rowstart[i] - 
            rowstart[i - 1]));
    return m;
  }

  index_t row_length(const index_t r) const {
    Assert( r < rows, "exceed range");
    return rowstart[r + 1] - rowstart[r];
  }

  std::size_t * get_rowstart_indices() {
    return rowstart;
  }

  const std::size_t * get_rowstart_indices() const {
    return rowstart;
  }

  index_t * get_column_numbers() {
    return colnums;
  }

  const index_t * get_column_numbers() const {
    return colnums;
  }

private:
  index_t max_dim; 
  index_t rows, cols; 
  index_t max_vec_len; 
  index_t max_row_length; 

  index_t * colnums; 
  std::size_t * rowstart;

  bool compressed; 

  bool store_diagonal_first_in_row; 
  template <typename, typename> friend class SparseMatrix;
};

template <typename ostream, typename Index>
ostream& operator << (ostream& os, const SparsityPatternBase<Index>& sp) {
  const auto rowstart = sp.get_rowstart_indices();
  const auto colnums = sp.get_column_numbers(); 

  if(rowstart == NULL || colnums == NULL) return os; 

  for(auto i = 0; i < sp.n_rows(); ++i) 
    for(auto j = rowstart[i]; j < rowstart[i + 1]; ++j) 
      if(colnums[j] != sp.invalid_entry) 
        os << "(" << i << ", " << colnums[j] << ")\n";
  return (os << std::flush);
}

namespace Utilities {
  template <typename Iterator, typename T, typename Comp>
    Iterator lower_bound(Iterator, Iterator, const T&, const Comp);

  template <typename Iterator, typename T>
    Iterator lower_bound(Iterator, Iterator, const T&);
};

#include "sparsity_pattern.templates.h"

class SparsityPattern : public SparsityPatternBase<unsigned int> {
  typedef SparsityPatternBase<unsigned int> base_t;
public:
  using base_t::base_t; 

  SparsityPattern() : base_t() {}

  SparsityPattern(const SparsityPattern& other) : base_t(other) {}

  SparsityPattern& operator = (const SparsityPattern& other) {
    this->operator= (other);
    return *this;
  }

  template <typename, typename> friend class SparseMatrix;
};

namespace Utilities {

  template <typename Iterator, typename T, typename Comp>
  inline Iterator
  lower_bound(Iterator first, Iterator last, 
      const T& val, const Comp comp){
    Assert(last - first >= 0, "the iterators do not exist");

    unsigned int len = static_cast<unsigned int>(last - first);
    if(len == 0) return first;

    while(true){
      if(len < 8){
        switch(len){
          case 7:
            if(!comp(*first, val))
              return first;
            ++first;
          case 6:
            if(!comp(*first, val))
              return first;
            ++first;
          case 5:
            if(!comp(*first, val))
              return first;
            ++first;
          case 4:
            if(!comp(*first, val))
              return first;
            ++first;
          case 3:
            if(!comp(*first, val))
              return first;
            ++first;
          case 2:
            if(!comp(*first, val))
              return first;
            ++first;
          case 1:
            if(!comp(*first, val))
              return first;
            return first + 1;
          default:
            std::abort();
        }
      }

      const unsigned int half = len >> 1;
      const Iterator middle = first + half;
      if(comp(*middle, val)){
        first = middle + 1;
        len -= half + 1;
      }
      else
        len = half;
    }
  }

  template <typename Iterator, typename T>
    inline Iterator
    lower_bound(Iterator first, Iterator last, const T& val){
      return Utilities::lower_bound(first, last, val, std::less<T>());
    }
}

#endif
