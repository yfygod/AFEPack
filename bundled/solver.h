/**
 * @file solver.h
 * @brief Some sparse solvers, maily borrow from dealii
 * @author Yang Fanyi <yfy__92@126.com>
 * @version 
 * @date 2019-11-04
 */

#ifndef __SPARSE_SOLVER_H_
#define __SPARSE_SOLVER_H_

#include <cmath>

#include "exception.h"
#include "vector.h"
#include "sparse_matrix.h"

class SolverControl;
template <typename VECTOR> class Solver;
template <typename VECTOR> class SolverCG;
template <typename VECTOR> class SolverBiCGSTAB;

class PreconditionIdentity;
template <typename Scalar, typename Index = unsigned int> class SparseLUDecomposition;
template <typename Scalar, typename Index = unsigned int> class SparseILU;

class SolverControl {
public:
  enum State {
    iterate = 0, 
    success = 1,
    failure = -1,
  };

  SolverControl(const int mi = 100, 
      const double t = 1.e-10,
      const bool mlh = false,
      const bool mlr = true) : 
    maxsteps(mi), tol(t), lvalue(1.e+100), lstep(0), 
    m_log_history(mlh), m_log_frequency(1), 
    m_log_result(mlr), history_data_enabled(false) {

    }

  virtual ~SolverControl() {}

  virtual State check(const unsigned int, const double);

  State last_check() const { return lcheck;}

  double initial_value() const { return initial_val;}

  double last_value() const { return lvalue; }

  unsigned int last_step() const { return lstep;}

  unsigned int max_steps() const { return maxsteps;}

  unsigned int set_max_steps(const unsigned int newval) {
    unsigned int old = maxsteps; 
    maxsteps = newval;
    return old;
  }

  double tolerance() const { return tol; }

  double set_tolerance(const double t) {
    double old = tol;
    tol = t; 
    return old;
  }

  void enable_history_data() { history_data_enabled = true;}

  double average_reduction() const {
    if(lstep == 0)
      return 0.;
    AssertExc(history_data_enabled, "require enable_history_data()");

    return std::pow(history_data[lstep]/history_data[0], 
        1./lstep);
  }

  double final_reduction() const { return step_reduction(lstep);}

  double step_reduction(const unsigned int step) const {
    AssertExc(history_data_enabled, "require  enable_history_data()");
    AssertExc(step <= lstep, "step is too much");

    return history_data[step]/history_data[step - 1];
  }

  void set_log_history(const bool);

  bool log_history() const { return m_log_history;}

  void log_history(const bool newval) { m_log_history = newval;}

  unsigned int log_frequency(unsigned int f) {
    if(f == 0)
      f = 1;

    unsigned int old = m_log_frequency;
    m_log_frequency = f;
    return old;
  }

  void set_log_result(const bool);

  bool log_result() const { return m_log_result;}

  void log_result(const bool newval) { m_log_result = newval;}

protected:
  unsigned int maxsteps;
  
  double tol;

  State lcheck;

  double initial_val;

  double lvalue;

  unsigned int lstep;

  bool m_log_history;

  unsigned int m_log_frequency;

  bool m_log_result;

  bool history_data_enabled;

  std::vector<double> history_data;
  
};


inline SolverControl::State
SolverControl::check(const unsigned int step,
    const double check_value) {

  if(step == 0) {
    initial_val = check_value;
    if(history_data_enabled)
      history_data.resize(maxsteps);
  }

  if(m_log_history && ((step % m_log_frequency) == 0))
    std::cerr << "Check: " << step << "\t" << 
      check_value << std::endl;

  lstep = step;
  lvalue = check_value;

  if(step == 0)
    if(m_log_result)
      std::cerr << "Starting value " << check_value << std::endl;

  if(history_data_enabled)
    history_data[step] = check_value;

  if(check_value <= tol) {
    if(m_log_result)
      std::cerr << "Convergence step " << step
        << "\tvalue " << check_value << std::endl;
    lcheck = success;
    return success;
  }

  if( (step >= maxsteps) || std::isnan(check_value)) {
    if(m_log_result)
      std::cerr << "Failure step " << step
        << "\tvalue " << check_value << std::endl;
    lcheck = failure;
    return failure;
  }

  lcheck = iterate;

  return iterate;
} 

template <typename VECTOR>
class Solver {
public:
  Solver(SolverControl& solver_control) : cntrl(solver_control) {}

  SolverControl& control() const { return cntrl;}

protected:
  SolverControl& cntrl;
};

#include "solver/precondition.h"
#include "solver/solver_cg.h"
#include "solver/solver_bicgstab.h"
#include "solver/solver_gmres.h"
#include "solver/sparse_decomposition.h"
#include "solver/sparse_ilu.h"

#endif
