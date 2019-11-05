template <typename Scalar, typename Index>
class SparseLUDecomposition : public SparseMatrix<Scalar, Index> {
protected:
  SparseLUDecomposition();

  SparseLUDecomposition(const SparsityPatternBase<Index>&);

public:
  typedef Index size_type;

  virtual ~SparseLUDecomposition() = 0;

  virtual void clear();

  struct AdditionalData {
  public:
    AdditionalData(const double = 0,
        const unsigned int = 0,
        const bool = false,
        const SparsityPatternBase<Index> *  = NULL);

    double strengthen_diagonal;

    unsigned int extra_off_diagonals;

    bool use_previous_sparsity;
    
    const SparsityPatternBase<Index> * use_this_sparsity;
  };

  template <typename Scalar2>
  void initialize (const SparseMatrix<Scalar2, Index> &matrix,
                   const AdditionalData parameters);

  void reinit(const SparsityPatternBase<Index>&);

  template <typename Scalar2>
  void decompose (const SparseMatrix<Scalar2, Index>&, 
      const double = 0.);

  virtual bool is_decomposed() const;

  bool empty() const;

protected:
  template <typename Scalar2>
  void copy_from(const SparseMatrix<Scalar2, Index>&);

  virtual void strengthen_diagonal_impl();

  virtual Scalar get_strengthen_diagonal(const Scalar, const size_type row) const;

  bool decomposed;

  double  strengthen_diagonal;

  std::vector<const size_type *> prebuilt_lower_bound;

private:

  void prebuild_lower_bound ();

  SparsityPatternBase<Index> *own_sparsity;

};

template <typename Scalar, typename Index>
inline Scalar SparseLUDecomposition<Scalar, Index>::
get_strengthen_diagonal(const Scalar ,
                        const size_type ) const
{
  return strengthen_diagonal;
}



template <typename Scalar, typename Index>
bool SparseLUDecomposition<Scalar, Index>::is_decomposed () const
{
  return decomposed;
}



template <typename Scalar, typename Index>
bool SparseLUDecomposition<Scalar, Index>::empty () const
{
  return SparseMatrix<Scalar, Index>::empty();
}



template <typename Scalar, typename Index>
SparseLUDecomposition<Scalar, Index>::AdditionalData::AdditionalData (
  const double strengthen_diag,
  const unsigned int extra_off_diag,
  const bool use_prev_sparsity,
  const SparsityPatternBase<Index> *use_this_spars):
  strengthen_diagonal(strengthen_diag),
  extra_off_diagonals(extra_off_diag),
  use_previous_sparsity(use_prev_sparsity),
  use_this_sparsity(use_this_spars)
{}

template <typename Scalar, typename Index>
SparseLUDecomposition<Scalar, Index>::SparseLUDecomposition() : 
  SparseMatrix<Scalar, Index>(),
  decomposed(false),
  own_sparsity(NULL){

}


template <typename Scalar, typename Index>
SparseLUDecomposition<Scalar, Index>::
SparseLUDecomposition (const SparsityPatternBase<Index> &sparsity) :
  SparseMatrix<Scalar, Index>(sparsity),
  decomposed(false),
  own_sparsity(0)
{}



template<typename Scalar, typename Index>
SparseLUDecomposition<Scalar, Index>::~SparseLUDecomposition()
{
  clear();
}

template<typename Scalar, typename Index>
void SparseLUDecomposition<Scalar, Index>::clear()
{
  decomposed = false;

  std::vector<const size_type *> tmp;
  tmp.swap(prebuilt_lower_bound);

  SparseMatrix<Scalar, Index>::clear();

  if (own_sparsity)
    {
      delete own_sparsity;
      own_sparsity = NULL;
    }
}



template<typename Scalar, typename Index>
template <typename Scalar2>
void SparseLUDecomposition<Scalar, Index>::initialize (
  const SparseMatrix<Scalar2, Index> &matrix,
  const AdditionalData data)
{
  const SparsityPatternBase<Index> &matrix_sparsity=matrix.get_sparsity_pattern();

  const SparsityPatternBase<Index> *sparsity_pattern_to_use = NULL;

  if (data.use_this_sparsity)
    sparsity_pattern_to_use = data.use_this_sparsity;
  else if (data.use_previous_sparsity &&
           !this->empty() &&
           (this->m()==matrix.m()))
    {
      // Use the sparsity that was
      // previously used. Scalarhis is
      // particularly useful when
      // matrix entries change but
      // not the sparsity, as for the
      // case of several Newton
      // iteration steps on an
      // unchanged grid.
      sparsity_pattern_to_use = &this->get_sparsity_pattern();
    }
  else if (data.extra_off_diagonals==0)
    {
      // Use same sparsity as matrix
      sparsity_pattern_to_use = &matrix_sparsity;
    }
  else
    {
      // Create new sparsity

      // for the case that
      // own_sparsity wasn't deleted
      // before (e.g. by clear()), do
      // it here
      if (own_sparsity)
        {
          // release the sparsity
          SparseMatrix<Scalar, Index>::clear();
          // delete it
          delete own_sparsity;
        }

      // and recreate
      own_sparsity = new SparsityPatternBase<Index>(matrix_sparsity,
                                         matrix_sparsity.max_entries_per_row()
                                         +2*data.extra_off_diagonals,
                                         data.extra_off_diagonals);
      own_sparsity->compress();
      sparsity_pattern_to_use = own_sparsity;
    }

  // now use this sparsity pattern
  Assert (sparsity_pattern_to_use->n_rows()==sparsity_pattern_to_use->n_cols(),
          ",,");
  decomposed = false;
  {
    std::vector<const size_type *> tmp;
    tmp.swap (prebuilt_lower_bound);
  }
  SparseMatrix<Scalar, Index>::reinit (*sparsity_pattern_to_use);
}


template<typename Scalar, typename Index>
template<typename Scalar2>
void
SparseLUDecomposition<Scalar, Index>::
decompose (const SparseMatrix<Scalar2, Index> &matrix,
           const double                    strengthen_diagonal)
{
  decomposed = false;

  this->strengthen_diagonal = strengthen_diagonal;
  prebuild_lower_bound ();
  copy_from (matrix);
  decomposed = true;
}

template <typename Scalar, typename Index>
void SparseLUDecomposition<Scalar, Index>::reinit (const SparsityPatternBase<Index> &sparsity)
{
  Assert (sparsity.n_rows() == sparsity.n_cols(), "error");
  decomposed = false;
  {
    std::vector<const size_type *> tmp;
    tmp.swap (prebuilt_lower_bound);
  }
  SparseMatrix<Scalar, Index>::reinit (sparsity);
}

template<typename Scalar, typename Index>
void
SparseLUDecomposition<Scalar, Index>::prebuild_lower_bound()
{
  const size_type *const
  column_Scalars = this->get_sparsity_pattern().get_column_numbers();
  const std::size_t *const
  rowstart_indices = this->get_sparsity_pattern().get_rowstart_indices();
  const size_type N = this->m();

  prebuilt_lower_bound.resize (N);

  for (size_type row=0; row<N; row++)
    {
      prebuilt_lower_bound[row]
        = Utilities::lower_bound (&column_Scalars[rowstart_indices[row]+1],
                                  &column_Scalars[rowstart_indices[row+1]],
                                  row);
    }
}



template <typename Scalar, typename Index>
template <typename Scalar2>
void
SparseLUDecomposition<Scalar, Index>::copy_from (const SparseMatrix<Scalar2, Index> &matrix)
{
  // check whether we use the same sparsity
  // pattern as the input matrix
  
  AssertCond(&this->get_sparsity_pattern() == &matrix.get_sparsity_pattern());

  if (&this->get_sparsity_pattern() == &matrix.get_sparsity_pattern())
    {
      const Scalar2 *input_ptr = matrix.get_nonzero_value();
      Scalar *this_ptr = this->get_nonzero_value();
      const Scalar *const end_ptr = this_ptr + this->n_nonzero_elements();
      if (std::is_same<Scalar2, Scalar>::value == true)
        std::memcpy (this_ptr, input_ptr, this->n_nonzero_elements()*sizeof(Scalar));
      else
        for ( ; this_ptr != end_ptr; ++input_ptr, ++this_ptr)
          *this_ptr = *input_ptr;
      return;
    }

  // preset the elements by zero. this needs to be written in a slightly
  // awkward way so that we find the corresponding function in the base class.
  SparseMatrix<Scalar, Index>::operator= (Scalar(0));

  // both allow more and less entries in the new matrix
/*  for (size_type row=0; row<this->m(); ++row)*/
    //{
      //typename SparseMatrix<Scalar, Index>::iterator index = this->begin(row);
      //typename SparseMatrix<Scalar2, Index>::const_iterator
      //in_index = matrix.begin(row);
      //index->value() = in_index->value();
      //++index, ++in_index;
      //while (index < this->end(row) && in_index < matrix.end(row))
        //{
          //while (index->column() < in_index->column() && index < this->end(row))
            //++index;
          //while (in_index->column() < index->column() && in_index < matrix.end(row))
            //++in_index;

          //index->value() = in_index->value();
          //++index, ++in_index;
        //}
    /*}*/
}



template <typename Scalar, typename Index>
void
SparseLUDecomposition<Scalar, Index>::strengthen_diagonal_impl ()
{
  for (size_type row=0; row<this->m(); ++row) {
    AssertCond (this->m() == this->n());
    //typename SparseMatrix<Scalar, Index>::iterator
    //diagonal_element = this->begin(row);

    const std::size_t * rowstart = this->get_rowstart_indices();
    const size_type * colnums = this->get_column_numbers();
    const Scalar* vals = this->get_nonzero_value();

    Scalar rowsum = 0;

    for(size_type j = rowstart[row] + 1; j < rowstart[row + 1]; ++j) {
      rowsum += std::fabs(vals[colnums[j]]);
    }

    /*      for (typename SparseMatrix<Scalar, Index>::iterator*/
    //p = diagonal_element + 1;
    //p != this->end(row); ++p)
    /*rowsum += std::fabs(p->value());*/

    this->diag_element(row) += this->get_strengthen_diagonal(rowsum, row)*rowsum;

    //diagonal_element->value() += this->get_strengthen_diagonal (rowsum, row)  *
    //rowsum;
  }
}


