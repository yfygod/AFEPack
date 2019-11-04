template <typename Scalar, typename Index>
SparseMatrix<Scalar, Index>& 
SparseMatrix<Scalar, Index>:: operator = (const Scalar d) {
  Assert(cols != NULL, "should be initialized first");
  Assert(cols->compressed || cols->empty(), "empty sparsity pattern");

  const auto nnzero = cols->n_nonzero_elements(); 

  std::fill(&val[0], &val[0] + nnzero, d);
  return *this;
}

template <typename Scalar, typename Index>
SparseMatrix<Scalar, Index>& 
SparseMatrix<Scalar, Index>:: operator += (const SparseMatrix& m) {
  AssertExc(cols == m.cols, "operator += must be used with the same SparsityPattern");
  AssertExc(val != m.val, "operator += must be used with different matrix");

  scalar_t * val_ptr = &val[0];
  scalar_t * end_ptr = &val[cols->n_nonzero_elements()];
  scalar_t * m_ptr = &m.val[0];

  for(; val_ptr != end_ptr; ++val_ptr, ++m_ptr)
    (*val_ptr) += (*m_ptr);
  return *this;
}

template <typename Scalar, typename Index>
void SparseMatrix<Scalar, Index>::reinit(
    const SparsityPatternBase<Index>& sp) {
  cols = &sp; 

  if(cols->empty()){
    if(val != NULL) delete[] val;
    val = NULL;
    max_len = 0;
    return;
  }

  const std::size_t N = cols->n_nonzero_elements();
  if(N > max_len || max_len == 0){
    if(val != NULL) delete[] val;
    val = new scalar_t[N];
    max_len = N;
  }

  std::memset(&val[0], 0, N*sizeof(scalar_t));
}

template <typename Scalar, typename Index>
void SparseMatrix<Scalar, Index>::clear() {
  cols = NULL;
  if(val) delete val;
  val = NULL;
  max_len = 0;
}

template <typename Scalar, typename Index>
Index SparseMatrix<Scalar, Index>::n_actually_nonzero_elements(
    const double threshold) const {
  Assert(threshold >= 0, "threshold should be positive");
  //Assert(cols != NULL, "should be initialized first");
  
  if(cols == NULL) return 0;
  const index_t nnz_alloc = n_nonzero_elements();
  index_t nnz = std::count_if(val, val + nnz_alloc, [threshold](const scalar_t n) { return n > threshold; });
  return nnz;
}

template <typename Scalar, typename Index>
void SparseMatrix<Scalar, Index>::set (const index_t i, 
    const index_t j, const scalar_t value) {
  const index_t index = cols->operator()(i, j);

  if (index == cols->invalid_entry)
  {
    AssertExc ((index != cols->invalid_entry) ||
        (value == 0.),
        "(i, j) should be added first");
    return;
  }

  val[index] = value;
}

template <typename Scalar, typename Index>
void SparseMatrix<Scalar, Index>::set_row(const index_t i, 
    const scalar_t value) {
  scalar_t * val_ptr = &val[cols->rowstart[i]];
  const scalar_t * const end_ptr =&val[cols->rowstart[i + 1]];

  while(val_ptr != end_ptr)
    *val_ptr ++ = value;
}

template <typename Scalar, typename Index>
void SparseMatrix<Scalar, Index>::add(const index_t i, 
    const index_t j, const scalar_t value) {
  if (value == 0)
    return;

  const index_t index = cols->operator()(i, j);

  if (index == cols->invalid_entry)
  {
    //AssertCond(false);
    AssertExp2 ((index != cols->invalid_entry), i, j);
    return;
  }

  val[index] += value;
}

template <typename Scalar, typename Index>
void SparseMatrix<Scalar, Index>::set_zero() {
  scalar_t value = 0.;
  (*this) = value;
}

template <typename Scalar, typename Index>
SparseMatrix<Scalar, Index>& 
SparseMatrix<Scalar, Index>::operator *= (const Scalar factor) {
  Assert(cols != NULL, "initialized first"); 
  Assert(val != NULL, "initialized first"); 
  scalar_t * val_ptr = &val[0];
  const scalar_t * const end_ptr = &val[cols->n_nonzero_elements()];
  while(val_ptr != end_ptr)
    *val_ptr++ *= factor;

  return *this;
}

template <typename Scalar, typename Index>
SparseMatrix<Scalar, Index>& 
SparseMatrix<Scalar, Index>::operator /= (const Scalar factor) {
  Assert(cols != NULL, "initialized first"); 
  Assert(val != NULL, "initialized first"); 

  const scalar_t factor_inv = 1./factor;

  scalar_t * val_ptr = &val[0];
  const scalar_t * const end_ptr = &val[cols->n_nonzero_elements()];
  while(val_ptr != end_ptr)
    *val_ptr++ *= factor_inv;

  return *this;
}

template <typename Scalar, typename Index>
template <typename OUTVEC, typename INVEC>
void SparseMatrix<Scalar, Index>::vmult(OUTVEC& out, 
    const INVEC& in) const {
  Assert(cols != NULL && val != NULL, "not initialized");
  Assert(&out != &in, "...");
  AssertExc(out.size() == this->m(), "the size of vector should be equal to row");
  AssertExc(in.size() == this->n(), "the size of vector should be equal to col");

  const std::size_t * rowstart = cols->rowstart;
  const index_t * colnums = cols->colnums;

  out.set_zero();

  for(index_t i = 0; i < cols->rows; ++i)
    for(std::size_t j = rowstart[i]; j < rowstart[i + 1]; ++j)
      out[i] += val[j]*in[colnums[j]];
}

template <typename Scalar, typename Index>
template <typename OUTVEC, typename INVEC>
void SparseMatrix<Scalar, Index>::Tvmult(OUTVEC& dst, 
    const INVEC& src) const {
  Assert(cols != NULL && val != NULL, "not initialized");
  Assert(&dst != &src, "...");
  AssertExc(dst.size() == this->n(), "the size of vector should be equal to col");
  AssertExc(src.size() == this->m(), "the size of vector should be equal to row");

  dst.set_zero();

  for (index_t i=0; i<m(); i++) {
    for (std::size_t j=cols->rowstart[i]; j<cols->rowstart[i+1] ; j++) {
      const index_t p = cols->colnums[j];
      dst(p) += val[j] * src(i);
    }
  }
}


template <typename Scalar, typename Index>
template <typename VEC> 
double SparseMatrix<Scalar, Index>::residual(const VEC& v1, 
    const VEC& v2) const {
  VEC t(v2.size()); vmult(t, v1); 
  t -= v2; 
  return t.template lpNorm<2>();
}

template <typename Scalar, typename Index>
template <typename INVEC1, typename INVEC2>
Scalar SparseMatrix<Scalar, Index>::vTmultv(const INVEC1& vec1, 
    const INVEC2& vec2) const {
  Assert(cols != NULL && val != NULL, "not initialized");
  AssertExc(vec1.size() == this->m(), "the size of vector should be equal to row");
  AssertExc(vec2.size() == this->n(), "the size of vector should be equal to col");

  const std::size_t * rowstart = cols->rowstart;
  const index_t * colnums = cols->colnums;

  scalar_t res(0.);

  for(index_t i = 0; i < cols->rows; ++i)
    for(std::size_t j = rowstart[i]; j < rowstart[i + 1]; ++j)
      res += vec1[i]*val[j]*vec2[colnums[j]];

  return res;
}

template <typename Scalar, typename Index>
void SparseMatrix<Scalar, Index>::print(const std::string& file, 
    int precision) const {
  std::ofstream fout(file); fout.precision(precision);
  fout << (*this);
  fout.close();
}
