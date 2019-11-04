#include <math.h>

void to3d_local_to_global(const double * lp, 
    const double ** lv, 
    const double ** gv, 
    double * gp) {
  double lambda[2];
  lambda[0] = (lv[1][0] - lp[0])/(lv[1][0] - lv[0][0]);
  lambda[1] = (lp[0] - lv[0][0])/(lv[1][0] - lv[0][0]);

  gp[0] = lambda[0]*gv[0][0] + lambda[1]*gv[1][0];
  gp[1] = lambda[0]*gv[0][1] + lambda[1]*gv[1][1];
  gp[2] = lambda[0]*gv[0][2] + lambda[1]*gv[1][2];
}

void to3d_global_to_local(const double * gp, 
    const double ** lv, 
    const double ** gv, 
    double * lp) {
  double lambda[2];
	lambda[1] = sqrt(
        (gv[1][0] - gv[0][0])*(gv[1][0] - gv[0][0]) 
			+ (gv[1][1] - gv[0][1])*(gv[1][1] - gv[0][1]) 
      + (gv[1][2] - gv[0][2])*(gv[1][2] - gv[0][2]));

	lambda[0] = sqrt(
        (gv[1][0] - gp[0])*(gv[1][0] - gp[0])
			+ (gv[1][1] - gp[1])*(gv[1][1] - gp[1])
      + (gv[1][2] - gp[2])*(gv[1][2] - gp[2])
      )/lambda[1];

  lambda[1] = 1. - lambda[0];
  lp[0] = lambda[0]*lv[0][0] + lambda[1]*lv[1][0];
}

double to3d_local_to_global_jacobian(const double * lp, 
		const double ** lv, 
		const double ** gv)
{
	double gl = sqrt(
        (gv[1][0] - gv[0][0])*(gv[1][0] - gv[0][0])
			+ (gv[1][1] - gv[0][1])*(gv[1][1] - gv[0][1])
      + (gv[1][2] - gv[0][2])*(gv[1][2] - gv[0][2])
      );
  return gl/fabs(lv[1][0] - lv[0][0]);
}

double to3d_global_to_local_jacobian(const double * gp, 
		const double ** lv,
		const double ** gv)
{
	double gl = sqrt(
        (gv[1][0] - gv[0][0])*(gv[1][0] - gv[0][0])
			+ (gv[1][1] - gv[0][1])*(gv[1][1] - gv[0][1])
      + (gv[1][2] - gv[0][2])*(gv[1][2] - gv[0][2])
      );
	return fabs(lv[1][0] - lv[0][0])/gl;
}
