namespace Impl {
  template <typename VECTOR>
  class TmpVectors {
  public:
    TmpVectors(const unsigned int max_size) : data(max_size, NULL), 
                                              offset(0) {}

    ~TmpVectors() {
      for(VECTOR * ptr : data) 
        if(ptr != NULL)
          delete ptr;
    }

    VECTOR& operator[] (const unsigned int i) const {
      VECTOR * ptr = data[i - offset];
      //AssertCond(ptr != NULL);
      return *ptr;
    }

    VECTOR& operator() (const unsigned int i, 
        const VECTOR& temp) {
      VECTOR * ptr = data[i - offset];
      if(ptr == NULL) 
        data[i - offset] = new VECTOR(temp);
      return *data[i - offset];
    }

  private:
    std::vector<VECTOR *> data;
    unsigned int offset;
  };
};

template <typename VECTOR>
class SolverGMRES : public Solver<VECTOR> {
public:
  struct AdditionalData {
    AdditionalData(const unsigned int m = 30, 
        const bool r = false, 
        const bool u = true,
        const bool f = false) : max_n_tmp_vectors(m), 
    right_preconditioning(r), 
    use_default_residual(u),
    force_re_orthogonalization(f) {}

    unsigned int max_n_tmp_vectors;

    bool right_preconditioning;

    bool use_default_residual;

    bool force_re_orthogonalization;
  };

  SolverGMRES(SolverControl& cn, 
      const AdditionalData& data = AdditionalData()) : 
    Solver<VECTOR>(cn), additional_data(data) {}

  template <typename MATRIX, typename PRE>
  void solve(const MATRIX&, VECTOR&, const VECTOR&, const PRE&);

private:
  virtual double criterion() { return 0.;}

  void givens_rotation (VECTOR&h,  VECTOR &b,
                        VECTOR&ci, VECTOR &si,
                        int col) const;

  static double
  modified_gram_schmidt (const Impl::TmpVectors<VECTOR> &orthogonal_vectors,
                         const unsigned int  dim,
                         const unsigned int  accumulated_iterations,
                         VECTOR             &vv,
                         VECTOR     &h,
                         bool               &re_orthogonalize);

private:
  AdditionalData additional_data;
  FullMatrix<double> H, H1;

};

template <typename VECTOR>
void SolverGMRES<VECTOR>::givens_rotation(
    VECTOR& h, VECTOR& b, VECTOR& ci, VECTOR& si, int col) const {
  for (int i=0 ; i<col ; i++) {
    const double s = si(i);
    const double c = ci(i);
    const double dummy = h(i);
    h(i)   =  c*dummy + s*h(i+1);
    h(i+1) = -s*dummy + c*h(i+1);
  };

  const double r = 1./std::sqrt(h(col)*h(col) + h(col+1)*h(col+1));
  si(col) = h(col+1) *r;
  ci(col) = h(col)   *r;
  h(col)  =  ci(col)*h(col) + si(col)*h(col+1);
  b(col+1)= -si(col)*b(col);
  b(col) *=  ci(col);
}

template <typename VECTOR>
double SolverGMRES<VECTOR>::modified_gram_schmidt(
    const Impl::TmpVectors<VECTOR>& orthogonal_vectors,
    const unsigned int dim, 
    const unsigned int  accumulated_iterations,
    VECTOR             &vv,
    VECTOR     &h,
    bool               &re_orthogonalize) {
  const unsigned int inner_iteration = dim - 1;
  double norm_vv_start = 0;

  if (re_orthogonalize == false && inner_iteration % 5 == 4)
    norm_vv_start = vv.l2_norm();

  for (unsigned int i=0 ; i<dim ; ++i) {
    h(i) = vv.dot(orthogonal_vectors[i]);
    vv.add(-h(i), orthogonal_vectors[i]);
  };

  if (re_orthogonalize == false && inner_iteration % 5 == 4) {
    const double norm_vv = vv.l2_norm();
    if (norm_vv > 10. * norm_vv_start *
        std::sqrt(std::numeric_limits<typename VECTOR::scalar_t>::epsilon()))
      return norm_vv;

    else {
      re_orthogonalize = true;
      std::cerr << "Re-orthogonalization enabled at step "
        << accumulated_iterations << std::endl;
    }
  }

  if (re_orthogonalize == true)
    for (unsigned int i=0 ; i<dim ; ++i) {
      double htmp = vv.dot(orthogonal_vectors[i]);
      h(i) += htmp;
      vv.add(-htmp, orthogonal_vectors[i]);
    }

  return vv.l2_norm();
}

template <typename VECTOR>
template <typename MATRIX, typename PRE>
void SolverGMRES<VECTOR>::solve(const MATRIX& A,
    VECTOR& x, const VECTOR& b, const PRE& precondition) {
  const unsigned int n_tmp_vectors = additional_data.max_n_tmp_vectors;
  Impl::TmpVectors<VECTOR> tmp_vectors(n_tmp_vectors);

  unsigned int accumulated_iterations = 0;
  H.reinit(n_tmp_vectors, n_tmp_vectors-1);

  VECTOR gamma(n_tmp_vectors),
        ci   (n_tmp_vectors-1),
        si   (n_tmp_vectors-1),
        h    (n_tmp_vectors-1);

  unsigned int dim = 0;

  SolverControl::State iteration_state = SolverControl::iterate;

  const bool left_precondition = !additional_data.right_preconditioning;

  const bool use_default_residual = additional_data.use_default_residual;

  VECTOR& v = tmp_vectors(0, x);
  VECTOR& p = tmp_vectors(n_tmp_vectors - 1, x);

/*  VECTOR *r      = NULL;*/
  //VECTOR *x_     = NULL;
  /*VECTOR *gamma_ = NULL;*/

  bool re_orthogonalize = additional_data.force_re_orthogonalization;

  do {
    h.reinit(n_tmp_vectors - 1);

    if (left_precondition) {
      A.vmult(p,x);
      p.sadd(-1.,1.,b);
      precondition.vmult(v,p);
    }
    else {
      A.vmult(v,x);
      v.sadd(-1.,1.,b);
    };

    double rho = v.l2_norm();

    iteration_state = this->control().check (accumulated_iterations, rho);

    if (iteration_state != SolverControl::iterate)
      break;


    //std::cerr << "default_res=" << rho << std::endl;

    /*    if (left_precondition) {*/
    //A.vmult(*r,x);
    //r->sadd(-1.,1.,b);
    //}
    //else
    //precondition.vmult(*r,v);

    //double res = r->l2_norm();
    //iteration_state = this->control().check (
    //accumulated_iterations, res);

    //if (iteration_state != SolverControl::iterate) {
    //if(r != NULL) delete r;
    //if(x_ != NULL) delete x_;
    //if(gamma_ != NULL) delete gamma_;
    //break;
    /*}*/

    gamma(0) = rho;

    v *= 1./rho;
    for (unsigned int inner_iteration=0;
        ((inner_iteration < n_tmp_vectors-2)
         &&
         (iteration_state==SolverControl::iterate));
        ++inner_iteration) {
      ++ accumulated_iterations;

      VECTOR &vv = tmp_vectors(inner_iteration+1, x);

      if (left_precondition) {
        A.vmult(p, tmp_vectors[inner_iteration]);
        precondition.vmult(vv,p);
      }
      else {
        precondition.vmult(p, tmp_vectors[inner_iteration]);
        A.vmult(vv,p);
      }; 

      dim = inner_iteration+1;
      const double s = modified_gram_schmidt(tmp_vectors, dim,
          accumulated_iterations,
          vv, h, re_orthogonalize);

      h(inner_iteration+1) = s;

      if (std::isfinite(1./s)) vv *= 1./s;

      givens_rotation(h,gamma,ci,si,inner_iteration);

      for (unsigned int i=0; i<dim; ++i) H(i,inner_iteration) = h(i);

      rho = std::fabs(gamma(dim));
      iteration_state = this->control().check (accumulated_iterations, rho);

    }

    h.reinit(dim);
    H1.reinit(dim+1,dim);

    for (unsigned int i=0; i<dim+1; ++i)
      for (unsigned int j=0; j<dim; ++j)
        H1(i,j) = H(i,j);

    H1.backward(h,gamma);

    if (left_precondition)
      for (unsigned int i=0 ; i<dim; ++i)
        x.add(h(i), tmp_vectors[i]);
    else {
      p = 0.;
      for (unsigned int i=0; i<dim; ++i)
        p.add(h(i), tmp_vectors[i]);
      precondition.vmult(v,p);
      x.add(1.,v);
    };
  } while(iteration_state == SolverControl::iterate);

  if (this->control().last_check() != SolverControl::success)
    std::cerr << "No convergence! " << "\n" << 
      "last step " << this->control().last_step() << "\n" <<
      "last check " << this->control().last_check() << std::endl;
}
