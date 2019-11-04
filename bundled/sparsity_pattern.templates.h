template <typename Index>
void SparsityPatternBase<Index>::reinit(const index_t m, 
    const index_t n, 
    const index_t max_per_row) {
  const std::vector<index_t> row_lengths(m, max_per_row);
  reinit(m, n, row_lengths);
}

template <typename Index>
void SparsityPatternBase<Index>::reinit(const index_t m, 
    const index_t n, 
    const std::vector<index_t>& row_lengths) {
  AssertExc(row_lengths.size() == m, "wrong size of row_lengths");

  rows = m;  cols = n;

  if((m == 0) || (n == 0)){
    if(rowstart) delete[] rowstart;
    if(colnums) delete[] colnums;

    rowstart = NULL;
    colnums = NULL;
    max_vec_len = max_dim = rows = cols = 0;
    max_row_length = 0;
    compressed = false;
    return;
  }

  store_diagonal_first_in_row = (m == n);

  std::size_t vec_len = 0;
  for (index_t i = 0; i < m; ++i)
    vec_len += std::min(static_cast<index_t>(store_diagonal_first_in_row ?
                                               std::max(row_lengths[i], 1U) :
                                               row_lengths[i]),
                        n);

  if (vec_len == 0)
  {
    vec_len = 1;
    if (colnums)
    {
      delete[] colnums;
      colnums = 0;
    }

    max_vec_len = vec_len;
    colnums = new index_t[max_vec_len];
  }

  max_row_length = (row_lengths.size() == 0 ?
      0 :
      std::min (static_cast<index_t>(*std::max_element(row_lengths.begin(),
            row_lengths.end())),
        n));

  if (store_diagonal_first_in_row && (max_row_length==0) && (m!=0))
    max_row_length = 1;

  if (rows > max_dim)
  {
    if (rowstart)
    {
      delete[] rowstart;
      rowstart = 0;
    }

    max_dim = rows;
    rowstart = new std::size_t[max_dim+1];
  }

  // allocate memory for the column
  // numbers if necessary
  if (vec_len > max_vec_len)
  {
    if (colnums)
    {
      delete[] colnums;
      colnums = 0;
    }

    max_vec_len = vec_len;
    colnums = new index_t[max_vec_len];
  }

  // set the rowstart array
  rowstart[0] = 0;
  for (index_t i = 1; i <= rows; ++i)
    rowstart[i] = rowstart[i - 1] +
      (store_diagonal_first_in_row ?
       std::max(std::min(static_cast<index_t>(row_lengths[i - 1]),n),
         static_cast<index_t> (1U)) :
       std::min(static_cast<index_t>(row_lengths[i - 1]),n));
  Assert((rowstart[rows]==vec_len)
      ||
      ((vec_len == 1) && (rowstart[rows] == 0)),
      "internal error");

  std::fill_n (&colnums[0], vec_len, invalid_entry);

  if (store_diagonal_first_in_row)
    for (index_t i = 0; i < rows; i++)
      colnums[rowstart[i]] = i;

  compressed = false;
}

template <typename Index>
void SparsityPatternBase<Index>::add(const index_t i, 
    const index_t j) {
  Assert((rowstart != NULL) && (colnums != NULL), "empty sparse matrix");
  Assert(compressed == false, "compressed matrix");
  AssertExc(i < rows, "range error");
  AssertExc(j < cols, "range error");

  std::size_t k = rowstart[i];

  for(; k < rowstart[i + 1]; ++k){
    if(colnums[k] == j) return;
    if(colnums[k] == invalid_entry){
      colnums[k] = j;
      return;
    }
  }

  AssertExc(k < rowstart[i + 1], "run out space!");
}

template <typename Index>
void SparsityPatternBase<Index>::compress() {
  Assert((rowstart!=0) && (colnums!=0), "empty sparse matrix");

  if (compressed)
    return;

  index_t next_free_entry = 0,
            next_row_start  = 0,
            row_length      = 0;


  const std::size_t nonzero_elements
    = std::count_if (&colnums[rowstart[0]],
        &colnums[rowstart[rows]],
        std::bind2nd(std::not_equal_to<index_t>(), invalid_entry));

  index_t *new_colnums = new index_t[nonzero_elements];

  std::vector<index_t> tmp_entries (max_row_length);

  for (index_t line=0; line<rows; ++line)
  {
    row_length = 0;
    for (index_t j=rowstart[line]; j<rowstart[line+1]; ++j,++row_length)
      if (colnums[j] != invalid_entry)
        tmp_entries[row_length] = colnums[j];
      else
        break;

    if (row_length > 1)
      std::sort ((store_diagonal_first_in_row)
          ? tmp_entries.begin()+1
          : tmp_entries.begin(),
          tmp_entries.begin()+row_length);

    for (index_t j=0; j<row_length; ++j)
      new_colnums[next_free_entry++] = tmp_entries[j];

    rowstart[line] = next_row_start;
    next_row_start = next_free_entry;

    Assert ((!store_diagonal_first_in_row) ||
        (new_colnums[rowstart[line]] == line),
        "internal error");

    Assert ((rowstart[line] == next_row_start)
        ||
        (std::find (&new_colnums[rowstart[line]+1],
                    &new_colnums[next_row_start],
                    new_colnums[rowstart[line]]) ==
         &new_colnums[next_row_start]),
        "internal error");
    Assert ((rowstart[line] == next_row_start)
        ||
        (std::adjacent_find(&new_colnums[rowstart[line]+1],
                            &new_colnums[next_row_start]) ==
         &new_colnums[next_row_start]),
        "internal error");
  };

  Assert (next_free_entry == nonzero_elements,
      "internal error");

  rowstart[rows] = next_row_start;

  delete[] colnums;
  colnums = new_colnums;

  max_vec_len = nonzero_elements;

  compressed = true;

}

template <typename Index>
Index SparsityPatternBase<Index>::operator() (const index_t i, 
    const index_t j) const {
  Assert((rowstart != NULL) && (colnums != NULL), "empty sparse matrix");
  Assert(i < rows, "i should be smaller than rows");
  Assert(j < cols, "j should be smaller than cols");
  Assert(compressed, "the operator() should be used after compressing");

  if (rowstart[i] == rowstart[i+1])
    return invalid_entry;

  if (store_diagonal_first_in_row && (i==j))
    return rowstart[i];

  const index_t *sorted_region_start = (store_diagonal_first_in_row ?
                                          &colnums[rowstart[i]+1] :
                                          &colnums[rowstart[i]]);
  const index_t *const p
    = Utilities::lower_bound<const index_t *> (sorted_region_start,
                                                 &colnums[rowstart[i+1]],
                                                 j);
  if ((p != &colnums[rowstart[i+1]])  &&  (*p == j))
    return (p - &colnums[0]);
  else
    return invalid_entry;

}
