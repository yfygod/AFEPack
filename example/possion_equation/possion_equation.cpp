////////////////////////////////////////////////////////////////////////////////////////////
// main1.cpp :
//

#include <iostream>
#include <fstream>
#include <sstream>

#include "Geometry.h"
#include "TemplateElement.h"
#include "FEMSpace.h"
#include "Operator.h"
#include "Functional.h"
#include "EasyMesh.h"
#include "AMGSolver.h"
#include "bundled/solver.h"

#define PI (4.0*atan(1.0))

double u(const double *);
double f(const double *);

using namespace AFEPack;

int main(int argc, char * argv[])
{

  HGeometryTree<2> htree;
  htree.readEasyMesh(argv[1]);
  IrregularMesh<2> ir_mesh;
  ir_mesh.reinit(htree);
  ir_mesh.globalRefine(atoi(argv[2]));
  ir_mesh.semiregularize();
  ir_mesh.regularize(true);

  Mesh<2,2> mesh = ir_mesh.regularMesh();

  std::stringstream tmp_dof, bas_fun;
  tmp_dof << "triangle." << argv[3] << ".tmp_dof";
  bas_fun << "triangle." << argv[3] << ".bas_fun";


  TemplateGeometry<2>	triangle_template_geometry;
  triangle_template_geometry.readData("triangle.tmp_geo");
  CoordTransform<2,2>	triangle_coord_transform;
  triangle_coord_transform.readData("triangle.crd_trs");
  TemplateDOF<2>	triangle_template_dof(triangle_template_geometry);
  triangle_template_dof.readData(tmp_dof.str());
  BasisFunctionAdmin<double,2,2> triangle_basis_function(triangle_template_dof);
  triangle_basis_function.readData(bas_fun.str());

  std::vector<TemplateElement<double,2,2> > template_element(1);
  template_element[0].reinit(triangle_template_geometry,
			     triangle_template_dof,
			     triangle_coord_transform,
			     triangle_basis_function);

  FEMSpace<double,2> fem_space(mesh, template_element);
	
  int n_element = mesh.n_geometry(2);
  fem_space.element().resize(n_element);
  for (int i = 0;i < n_element;i ++)
    fem_space.element(i).reinit(fem_space,i,0);

  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  StiffMatrix<2,double> stiff_matrix(fem_space);
  stiff_matrix.algebricAccuracy() = 10;
  stiff_matrix.build();

  FEMFunction<double,2> solution(fem_space);
  Vector<double> right_hand_side;
  Operator::L2Discretize(&f, fem_space, right_hand_side, 10);

  BoundaryFunction<double,2> boundary(BoundaryConditionInfo::DIRICHLET, 1, &u);
  BoundaryConditionAdmin<double,2> boundary_admin(fem_space);
  boundary_admin.add(boundary);
  boundary_admin.apply(stiff_matrix, solution, right_hand_side);

  //AMGSolver solver(stiff_matrix,1);
  //solver.solve(solution, right_hand_side, 1.0e-08, 200);	

  SolverControl cn; cn.log_history(true);
  cn.set_max_steps(1000);
  SolverBiCGSTAB<Vector<double>> solver(cn);
  SparseILU<double> pre; pre.initialize(stiff_matrix);
  solver.solve(stiff_matrix, solution, right_hand_side, pre);

  solution.writeOpenDXData("u.dx");
  double error = Functional::L2Error(solution, FunctionFunction<double>(&u), 10);
  //double error2 = Functional::H1SemiError(solution, FunctionFunction<double>(&u),10);
  std::cerr << "\nL2 error = " << error << std::endl;
  //std::cerr << "\nH1 error = " << error2 << std::endl;
  return 0;
};

double u(const double * p)
{
  return sin(1*PI*p[0]) * sin(2*PI*p[1]);
};

double f(const double * p)
{
  return 5*PI*PI*u(p);
};

//
// end of file
////////////////////////////////////////////////////////////////////////////////////////////
