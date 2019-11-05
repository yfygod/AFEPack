template <typename Scalar, typename Index>
class SparseILU : public SparseLUDecomposition<Scalar, Index> {
public:
  typedef Index size_type;
  
  SparseILU();
  
  SparseILU(const SparsityPatternBase<Index>&);

  typedef typename SparseLUDecomposition<Scalar, Index>::
    AdditionalData AdditionalData;

  template <typename Scalar2>
  void initialize(const SparseMatrix<Scalar2, Index>&,
      const AdditionalData& = AdditionalData());

  template <typename Scalar2>
  void decompose(const SparseMatrix<Scalar2, Index>&,
      const double  = 0.);

  template <typename Scalar2>
  void apply_decomposition(Vector<Scalar2>&,
      const Vector<Scalar2>&) const;

  template <typename Scalar2>
  void vmult(Vector<Scalar2>&, 
      const Vector<Scalar2>&) const;


  template <typename Scalar2>
  void Tvmult(Vector<Scalar2>&, 
      const Vector<Scalar2>&) const;

};

template <typename Scalar, typename Index>
template <typename Scalar2>
inline
void
SparseILU<Scalar, Index>::apply_decomposition (Vector<Scalar2>  &dst,
                              const Vector<Scalar2> &src) const
{
  vmult(dst, src);
}

template <typename Scalar, typename Index>
SparseILU<Scalar, Index>::SparseILU() {}

template <typename Scalar, typename Index>
SparseILU<Scalar, Index>::SparseILU(const SparsityPatternBase<Index>& sp) {
  SparseMatrix<Scalar, Index>::reinit(sp);
}

template <typename Scalar, typename Index>
template <typename Scalar2>
void 
SparseILU<Scalar, Index>::initialize(const SparseMatrix<Scalar2, Index>& matrix,
    const AdditionalData& data) {
  SparseLUDecomposition<Scalar, Index>::initialize(matrix, data);
  decompose(matrix, data.strengthen_diagonal);
}

template <typename Scalar, typename Index>
template <typename Scalar2>
void SparseILU<Scalar, Index>::decompose (const SparseMatrix<Scalar2, Index> &matrix,
                                   const double strengthen_diagonal)
{
  AssertExc (matrix.m()==matrix.n(), "error");
  AssertExc (this->m()==this->n(),   "error");
  AssertExc (matrix.m()==this->m(),  "error");

  AssertExc (strengthen_diagonal >= 0, "error");

  SparseLUDecomposition<Scalar, Index>::decompose (matrix, strengthen_diagonal);

  if(strengthen_diagonal > 0)
    this->strengthen_diagonal_impl();

  const SparsityPatternBase<Index>& sp = this->get_sparsity_pattern();
  const std::size_t * const ia = sp.get_rowstart_indices();
  const size_type * const ja = sp.get_column_numbers();

  Scalar* luval = this->SparseMatrix<Scalar, Index>::get_nonzero_value();
  const size_type N = this->m();
  size_type jrow = 0;

  std::vector<size_type> iw(N, static_cast<size_type>(-1));

  for(size_type k = 0; k < N; ++k) {
    const size_type j1 = ia[k], j2 = ia[k + 1] - 1;

    for(size_type j = j1; j <= j2; ++j)
      iw[ja[j]] = j;

    size_type j = j1 + 1;

    if(j > j2) 
      goto label_200;

label_150:

    jrow = ja[j];
    if(jrow >= k)
      goto label_200;

    {
      Scalar t1 = luval[j]*luval[ia[jrow]];
      luval[j] = t1;

      size_type jj = ia[jrow] + 1;
      while(ja[jj] < jrow)
        ++ jj;

      for(; jj < ia[jrow + 1]; ++jj) {
        const size_type jw = iw[ja[jj]];
        if(jw != static_cast<size_type>(-1))
          luval[jw] -= t1*luval[jj];
      }

      ++j;

      if(j <= j2)
        goto label_150;
    }

label_200:
    Assert( (jrow > k) || (j == ia[k + 1]), "error");
    Assert(luval[ia[k]] != 0, "error");

    luval[ia[k]] = 1./luval[ia[k]];

    for(size_type j = j1; j <= j2; ++j)
      iw[ja[j]] = static_cast<size_type>(-1);
  }
}

template <typename Scalar, typename Index>
template <typename Scalar2>
void 
SparseILU<Scalar, Index>::vmult(Vector<Scalar2>& dst,
    const Vector<Scalar2>& src) const {
  AssertCond(dst.size() == src.size());
  AssertCond(dst.size() == this->m());

  const size_type N = dst.size();
  const std::size_t* const rowstart_indices = 
    this->get_rowstart_indices();
  const size_type * const column_numbers = 
    this->get_column_numbers();

  dst = src;

  for(size_type row = 0; row < N; ++row) {
    const size_type * const rowstart = 
      &column_numbers[rowstart_indices[row] + 1];
    
    const size_type * const first_after_diagonal = 
      this->prebuilt_lower_bound[row];
    Scalar2 dst_row = dst(row);

    const Scalar* luval = this->get_nonzero_value() + 
      (rowstart - column_numbers);

    for(const size_type * col = rowstart; 
        col != first_after_diagonal; ++ col, ++ luval) 
      dst_row -= *luval*dst(*col);

    dst(row) = dst_row;
  }

  for(int row = N - 1; row >= 0; --row) {
    const size_type * const rowend = &column_numbers[
      rowstart_indices[row + 1]];
    const size_type * const first_after_diagonal = 
      this->prebuilt_lower_bound[row];

    Scalar2 dst_row = dst(row);

    const Scalar * luval = this->get_nonzero_value() + 
      (first_after_diagonal - column_numbers);

    for(const size_type * col = first_after_diagonal; 
        col != rowend; ++col, ++luval) 
      dst_row -= *luval*dst(*col);

    dst(row) = dst_row*this->diag_element(row);
  }
}

template <typename Scalar, typename Index>
template <typename Scalar2>
void
SparseILU<Scalar, Index>::Tvmult(Vector<Scalar2>& dst, 
    const Vector<Scalar2>& src) const {
  AssertCond(dst.size() == src.size());
  AssertCond(dst.size() == this->m());

  const size_type N = dst.size();
  const std::size_t* const rowstart_indices = 
    this->get_rowstart_indices();
  const size_type * const column_numbers = 
    this->get_sparsity_pattern();

  Vector<Scalar2> tmp(N);

  dst = src;
  
  for(size_type row = 0; row < N; ++row) {
    dst(row) -= tmp(row);
    
    dst(row) *= this->diag_element(row);
    const size_type * const rowend = &column_numbers[
      rowstart_indices[row + 1]];

    const size_type * const first_after_diagonal = 
      this->prebuilt_lower_bound[row];

    const Scalar2 dst_row = dst(row); 
    const Scalar * luval = this->get_nonzero_value() + 
      (first_after_diagonal - column_numbers);

    for(const size_type * col = first_after_diagonal; 
        col != rowend; ++col, ++luval)
      tmp(*col) += *luval*dst_row;
  }

  tmp = 0;
  
  for(int row = N - 1; row >= 0; -- row) {
    dst(row) -= tmp(row);

    const size_type * const rowstart = &column_numbers[
      rowstart_indices[row] + 1];

    const size_type * const first_after_diagonal = 
      this->prebuilt_lower_bound[row];
    
    const Scalar2 dst_row = dst(row);
    const Scalar* luval = this->get_nonzero_value() + 
      (rowstart - column_numbers);

    for(const size_type * col = rowstart;
        col != first_after_diagonal; ++col, ++ luval)
      tmp(*col) += *luval*dst_row;
  }
}
