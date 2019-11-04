///////////////////////////////////////////////////////////////////////////
// GmshMesh.h : by R.Lie
//

#ifndef _GmshMesh_h_
#define _GmshMesh_h_

#include <iostream>
#include <string>
#include <fstream>

#include "Geometry.h"

AFEPACK_OPEN_NAMESPACE

/**
 * This class provides facilities to asscess the mesh data file generated
 * by the mesh generator \p{gmsh}. For 3-dimensional only. Though we can
 * read in very flexible data format, we currently only use it to read in
 * pure tetrahedron mesh.
 */
class GmshMesh2D : public SimplestMesh<2,2>
{
public:
	typedef int GeometryType;
	static const GeometryType POINT			 = 15;
	static const GeometryType LINE			 = 1;
	static const GeometryType TRIANGLE		 = 2;
    static const GeometryType QUADRANGLE     = 3;
private:
	std::list<GeometryBM> node;
	std::list<GeometryBM> line;
	std::list<GeometryBM> surface;
public:
  GmshMesh2D();
  virtual ~GmshMesh2D();

  SimplestMesh<2, 2>& get_simplest_mesh() {
    return *this;
  }

  const SimplestMesh<2, 2>& get_simplest_mesh() const {
    return *this;
  }

  std::list<GeometryBM>& get_node_geometry_list() { 
    return node;
  }
  const std::list<GeometryBM>& get_node_geometry_list() const { 
    return node;
  }

  std::list<GeometryBM>& get_line_geometry_list() { 
    return line;
  }
  const std::list<GeometryBM>& get_line_geometry_list() const { 
    return line;
  }

  std::list<GeometryBM>& get_surface_geometry_list() { 
    return surface;
  }
  const std::list<GeometryBM>& get_surface_geometry_list() const { 
    return surface;
  }

public:
	void readData(const std::string&);
	virtual void generateMesh(Mesh<2,2>& m) override;
private:
	void base_generateMesh(Mesh<2,2>& m);
};


class GmshMesh3D : public SimplestMesh<3,3>
{
public:
	typedef int GeometryType;
	static const GeometryType POINT			= 15;
	static const GeometryType LINE			= 1;
	static const GeometryType TRIANGLE		= 2;
	static const GeometryType QUADRANGLE	= 3;
	static const GeometryType TETRAHEDRON	= 4;
	static const GeometryType HEXAHEDRON	= 5;
	static const GeometryType PRISM			= 6;
	static const GeometryType PYRAMID		= 7;
private:
	std::list<GeometryBM> node;
	std::list<GeometryBM> line;
	std::list<GeometryBM> surface;
public:
    GmshMesh3D();
    virtual ~GmshMesh3D();
public:
	void readData(const std::string&);
	virtual void generateMesh(Mesh<3,3>& m);

private:
	void base_generateMesh(Mesh<3,3>& m);
};




class GmshMesh : public SimplestMesh<3,3>
{
public:
	typedef int GeometryType;
	static const GeometryType POINT			= 15;
	static const GeometryType LINE			= 1;
	static const GeometryType TRIANGLE		= 2;
	static const GeometryType QUADRANGLE		= 3;
	static const GeometryType TETRAHEDRON		= 4;
	static const GeometryType HEXAHEDRON		= 5;
	static const GeometryType PRISM			= 6;
	static const GeometryType PYRAMID		= 7;
private:
	std::list<GeometryBM> node;
	std::list<GeometryBM> line;
  std::list<GeometryBM> surface;
public:
  GmshMesh();
  virtual ~GmshMesh();


public:
  void readData(const std::string&);
  virtual void generateMesh(Mesh<3,3>& m);
};

AFEPACK_CLOSE_NAMESPACE

#endif // _GmshMesh_h_

//
// end of file
///////////////////////////////////////////////////////////////////////////
