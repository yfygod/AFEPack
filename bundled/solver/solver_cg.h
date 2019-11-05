template <typename VECTOR>
class SolverCG : public Solver<VECTOR> {
  typedef Solver<VECTOR> base;
public:
  //typedef typename VECTOR::index_t size_type;

  SolverCG(SolverControl&);

  SolverCG(SolverControl&, const bool);

  virtual ~SolverCG();

  template <typename MATRIX, typename PRE> void
  solve(const MATRIX&, VECTOR&, const VECTOR&, 
      const PRE&);

  template <typename MATRIX, typename PRE> void
  solve(const MATRIX&, VECTOR&, const VECTOR&, 
      PRE&);

  void set_exact_residual(const bool);
  bool exact_residual() const;

protected:
  template <typename MATRIX>
  double criterion(const MATRIX&, const VECTOR&, const VECTOR&);

  virtual void print_vectors(const unsigned int,
      const VECTOR&, const VECTOR&, const VECTOR&) const;

  VECTOR Vr;
  VECTOR Vp;
  VECTOR Vz;
  VECTOR Vt;

  double res2;

  bool exact_res;

private:
  void cleanup();
};


template <typename VECTOR>
SolverCG<VECTOR>::SolverCG(SolverControl& cn) : 
  Solver<VECTOR>(cn), exact_res(false) {
}

template <typename VECTOR>
SolverCG<VECTOR>::SolverCG(SolverControl& cn, 
    const bool exact) : Solver<VECTOR>(cn), 
  exact_res(exact) {

}


template <typename VECTOR>
SolverCG<VECTOR>::~SolverCG() {}

template <typename VECTOR>
void
SolverCG<VECTOR>::set_exact_residual(const bool exact) {
  exact_res = exact;
}

template <typename VECTOR>
bool
SolverCG<VECTOR>::exact_residual() const {
  return exact_res;
}

template <typename VECTOR>
template <typename MATRIX>
double 
SolverCG<VECTOR>::criterion(const MATRIX& A, 
    const VECTOR& x,
    const VECTOR& b) {
  A.vmult(Vt, x);
  Vt -= b;
  return Vt.l2_norm();
}

template <typename VECTOR>
void
SolverCG<VECTOR>::cleanup() {

}

template <typename VECTOR>
void 
SolverCG<VECTOR>::print_vectors(const unsigned int, 
    const VECTOR&,
    const VECTOR&, 
    const VECTOR&) const {

}

template <typename VECTOR>
template <typename MATRIX, typename PRE>
void 
SolverCG<VECTOR>::solve(const MATRIX& A, 
    VECTOR& x, 
    const VECTOR& b,
    const PRE& precondition) {
  SolverControl::State conv = SolverControl::iterate;
  VECTOR& g = Vr; VECTOR& d = Vz; VECTOR& h = Vp;
  Vt.resize(x.size());
  g.resize(x.size()); d.resize(x.size()); h.resize(x.size());

  int it = 0;
  double res, gh, alpha, beta;
  
  if(!x.all_zero()) {
    A.vmult(g, x);
    g -= b;
  } else {
    g = b; 
    g *= -1.;
  }

  res = g.l2_norm();
  conv = this->control().check(0, res);
  if(conv) return;

  if(std::is_same<PRE, PreconditionIdentity>::value == false) {
    precondition.vmult(h, g);
    d = h; d *= -1;
    gh = g.dot(h);
  } else {
    d = g; d *= -1;
    gh = res*res;
  }

  while(conv == SolverControl::iterate) {
    ++ it;
    A.vmult(h, d);
    alpha = d.dot(h);

    AssertCond(alpha != 0);
    alpha = gh/alpha;

    g.add(alpha, h);
    x.add(alpha, d);
    res = g.l2_norm();

    print_vectors(it, x, g, d);

    if(exact_res)
      conv = this->control().check(it, criterion(A, x, b));
    else
      conv = this->control().check(it, res);
    if(conv != SolverControl::iterate)
      break;

    if(std::is_same<PRE, PreconditionIdentity>::value == false) {
      precondition.vmult(h, g);
      beta = gh;
      AssertCond(beta != 0.);

      gh = g.dot(h);
      beta = gh/beta;
      
      d *= beta; d -= h;
    } else {
      beta = gh;
      gh = res*res;
      beta = gh/beta;
      d *= beta; d -= g;
    }
  }

  if(this->control().last_check() != SolverControl::success) {
    std::cerr << "No convergence !" << "\n" <<
      "last step " << this->control().last_step() << "\n" << 
      "last value " << this->control().last_value() << std::endl;
  }
}

template <typename VECTOR>
template <typename MATRIX, typename PRE>
void 
SolverCG<VECTOR>::solve(const MATRIX& A, 
    VECTOR& x, 
    const VECTOR& b,
    PRE& precondition) {
  SolverControl::State conv = SolverControl::iterate;
  VECTOR& g = Vr; VECTOR& d = Vz; VECTOR& h = Vp;
  Vt.resize(x.size());
  g.resize(x.size()); d.resize(x.size()); h.resize(x.size());

  int it = 0;
  double res, gh, alpha, beta;
  
  if(!x.all_zero()) {
    A.vmult(g, x);
    g -= b;
  } else {
    g = b; 
    g *= -1.;
  }

  res = g.l2_norm();
  conv = this->control().check(0, res);
  if(conv) return;

  if(std::is_same<PRE, PreconditionIdentity>::value == false) {
    precondition.vmult(h, g);
    d = h; d *= -1;
    gh = g.dot(h);
  } else {
    d = g; d *= -1;
    gh = res*res;
  }

  while(conv == SolverControl::iterate) {
    ++ it;
    A.vmult(h, d);
    alpha = d.dot(h);

    AssertCond(alpha != 0);
    alpha = gh/alpha;

    g.add(alpha, h);
    x.add(alpha, d);
    res = g.l2_norm();

    print_vectors(it, x, g, d);

    if(exact_res)
      conv = this->control().check(it, criterion(A, x, b));
    else
      conv = this->control().check(it, res);
    if(conv != SolverControl::iterate)
      break;

    if(std::is_same<PRE, PreconditionIdentity>::value == false) {
      precondition.vmult(h, g);
      beta = gh;
      AssertCond(beta != 0.);

      gh = g.dot(h);
      beta = gh/beta;
      
      d *= beta; d -= h;
    } else {
      beta = gh;
      gh = res*res;
      beta = gh/beta;
      d *= beta; d -= g;
    }
  }

  if(this->control().last_check() != SolverControl::success) {
    std::cerr << "No convergence !" << "\n" <<
      "last step " << this->control().last_step() << "\n" << 
      "last value " << this->control().last_value() << std::endl;
  }
}
