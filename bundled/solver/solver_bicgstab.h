template <typename VECTOR>
class SolverBiCGSTAB : public Solver<VECTOR> {
public:
  SolverBiCGSTAB(SolverControl&);
  
  SolverBiCGSTAB(SolverControl&, const bool);

  virtual ~SolverBiCGSTAB();

  template <typename MATRIX, typename PRE> void
  solve(const MATRIX&,
      VECTOR&,
      const VECTOR&,
      PRE&);

  template <typename MATRIX, typename PRE> void
  solve(const MATRIX&,
      VECTOR&,
      const VECTOR&,
      const PRE&);

protected:
  template <typename MATRIX>
  double criterion(const MATRIX&, const VECTOR&, const VECTOR&);

  virtual void print_vectors(const unsigned int,
      const VECTOR&,
      const VECTOR&,
      const VECTOR&) const;

  VECTOR Vr, Vrbar, Vp, Vy, Vz, Vt, Vv;
  VECTOR * Vx;
  const VECTOR * Vb;

  double alpha, beta, omega, rho, rhobar, res;

  unsigned int step;

  bool exact_residual;

private:
  template <typename MATRIX>
  SolverControl::State start(const MATRIX& A);

  template <typename MATRIX, typename PRE> bool
  iterate(const MATRIX&, PRE&);

  template <typename MATRIX, typename PRE> bool
  iterate(const MATRIX&, const PRE&);
};

template <typename VECTOR>
SolverBiCGSTAB<VECTOR>::SolverBiCGSTAB(SolverControl& cn) : 
  Solver<VECTOR>(cn) , exact_residual(true){
}

template <typename VECTOR>
SolverBiCGSTAB<VECTOR>::SolverBiCGSTAB(SolverControl& cn, 
    const bool exact) : Solver<VECTOR>(cn), exact_residual(exact) {

}


template <typename VECTOR>
SolverBiCGSTAB<VECTOR>::~SolverBiCGSTAB() {}

template <typename VECTOR>
template <typename MATRIX>
double
SolverBiCGSTAB<VECTOR>::criterion(const MATRIX& A, 
    const VECTOR& x, 
    const VECTOR& b) {
  A.vmult(Vt, x);
  Vt -= b;
  res = Vt.l2_norm();

  return res;
}


template <typename VECTOR>
template <typename MATRIX>
SolverControl::State
SolverBiCGSTAB<VECTOR>::start(const MATRIX& A) {
  A.vmult(Vr, *Vx);
  Vr *= -1.; Vr += *Vb;
  res = Vr.l2_norm();

  return this->control().check(step, res);
}

template <typename VECTOR>
void
SolverBiCGSTAB<VECTOR>::print_vectors(const unsigned int,
    const VECTOR&,
    const VECTOR&,
    const VECTOR& ) const {

}

template <typename VECTOR>
template <typename MATRIX, typename PRE>
bool
SolverBiCGSTAB<VECTOR>::iterate(const MATRIX& A, 
    PRE& precondition) {
  SolverControl::State state = SolverControl::iterate;
  alpha = omega = rho = 1.;

  VECTOR& r = Vr; VECTOR& rbar = Vrbar; VECTOR& p = Vp;
  VECTOR& y = Vy; VECTOR& z = Vz; VECTOR& t = Vt;
  VECTOR& v = Vv;

  rbar = r;
  bool startup = true;

  do {
    ++ step;
    rhobar = r.dot(rbar);
    beta   = rhobar * alpha / (rho * omega);
    rho    = rhobar;
    if(startup == true) {
      p = r;
      startup = false; 
    } else {
      p *= beta; p += r; p -= (beta*omega*v);
    }

    precondition.vmult(y,p);
    A.vmult(v,y);
    rhobar = rbar.dot(v);

    alpha = rho/rhobar;

    if(std::fabs(alpha) > 1.e+10)
      return true;

    r -= (alpha*v);

    if(this->control().check(step, r.l2_norm()) == SolverControl::success) {
      Vx->add(alpha, y);
      print_vectors(step, *Vx, r, y);
      return false;
    }

    precondition.vmult(z, r);
    A.vmult(t, z);
    rhobar = t.dot(r);
    omega = rhobar/(t.dot(t));
    Vx->add(alpha, y); Vx->add(omega, z);
    r -= (omega*t);

    if(exact_residual)
      res = criterion(A, *Vx, *Vb);
    else
      res = r.l2_norm();

    state = this->control().check(step, res);
    print_vectors(step, *Vx, r, y);

  } while (state == SolverControl::iterate);

  return false;
}

template <typename VECTOR>
template <typename MATRIX, typename PRE>
void
SolverBiCGSTAB<VECTOR>::solve(const MATRIX& A,
    VECTOR& x,
    const VECTOR& b, 
    PRE& precondition) {
#define init(_r) _r.resize(x.size()); 
  init(Vr); init(Vrbar); init(Vp); init(Vy); init(Vz);
  init(Vt); init(Vv); 
#undef init
  Vx = &x; Vb = &b;

  step = 0;

  bool state;

  do {
    if(step != 0)
      std::cerr << "Restart step " << step << std::endl;
    if(start(A) == SolverControl::success)
      break;
    state = iterate(A, precondition);
  } while(state);

  if(this->control().last_check() != SolverControl::success)
    std::cerr << "No convergence! " << "\n" << 
      "last step " << this->control().last_step() << "\n" <<
      "last check " << this->control().last_check() << std::endl;
}


template <typename VECTOR>
template <typename MATRIX, typename PRE>
bool
SolverBiCGSTAB<VECTOR>::iterate(const MATRIX& A, 
    const PRE& precondition) {
  SolverControl::State state = SolverControl::iterate;
  alpha = omega = rho = 1.;

  VECTOR& r = Vr; VECTOR& rbar = Vrbar; VECTOR& p = Vp;
  VECTOR& y = Vy; VECTOR& z = Vz; VECTOR& t = Vt;
  VECTOR& v = Vv;

  rbar = r;
  bool startup = true;

  do {
    ++ step;
    rhobar = r.dot(rbar);
    beta   = rhobar * alpha / (rho * omega);
    rho    = rhobar;
    if(startup == true) {
      p = r;
      startup = false; 
    } else {
      p *= beta; p += r; p -= (beta*omega*v);
    }

    precondition.vmult(y,p);
    A.vmult(v,y);
    rhobar = rbar.dot(v);

    alpha = rho/rhobar;

    if(std::fabs(alpha) > 1.e+10)
      return true;

    r -= (alpha*v);

    if(this->control().check(step, r.l2_norm()) == SolverControl::success) {
      Vx->add(alpha, y);
      print_vectors(step, *Vx, r, y);
      return false;
    }

    precondition.vmult(z, r);
    A.vmult(t, z);
    rhobar = t.dot(r);
    omega = rhobar/(t.dot(t));
    Vx->add(alpha, y); Vx->add(omega, z);
    r -= (omega*t);

    if(exact_residual)
      res = criterion(A, *Vx, *Vb);
    else
      res = r.l2_norm();

    state = this->control().check(step, res);
    print_vectors(step, *Vx, r, y);

  } while (state == SolverControl::iterate);

  return false;
}

template <typename VECTOR>
template <typename MATRIX, typename PRE>
void
SolverBiCGSTAB<VECTOR>::solve(const MATRIX& A,
    VECTOR& x,
    const VECTOR& b, 
    const PRE& precondition) {
#define init(_r) _r.resize(x.size()); 
  init(Vr); init(Vrbar); init(Vp); init(Vy); init(Vz);
  init(Vt); init(Vv); 
#undef init
  Vx = &x; Vb = &b;

  step = 0;

  bool state;

  do {
    if(step != 0)
      std::cerr << "Restart step " << step << std::endl;
    if(start(A) == SolverControl::success)
      break;
    state = iterate(A, precondition);
  } while(state);

  if(this->control().last_check() != SolverControl::success)
    std::cerr << "No convergence! " << "\n" << 
      "last step " << this->control().last_step() << "\n" <<
      "last check " << this->control().last_check() << std::endl;
}

