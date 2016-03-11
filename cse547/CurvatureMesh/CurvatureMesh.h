/*! \file GraphViewerDynamicMesh.h
*   \brief Graph Viewer DynamicMesh
*   \author David Gu
*   \date   documented on 01/10/2014
*
*   Dynamic Mesh for viewer 
*/
#ifndef  _CURVATURE_MESH_H_
#define  _CURVATURE_MESH_H_

#include <map>
#include <vector>

#include "Mesh/Vertex.h"
#include "Mesh/Edge.h"
#include "Mesh/Face.h"
#include "Mesh/HalfEdge.h"

#include "Mesh/BaseMesh.h"
#include "Mesh/boundary.h"
#include "Mesh/iterators.h"
#include "Parser/parser.h"

#ifndef PI
#define PI 3.1415926535
#endif

namespace MeshLib
{
/*! \brief CCurvatureFace class
*
*	Face class for viewer
*   Trait: face normal
*/
  class CCurvatureFace: public CFace
  {
  protected:
    /*! face normal */
    CPoint m_normal;
	/*! face area   */
	double m_area;

  public:
    /*! face normal */
	  CPoint & normal() { return m_normal; };

	  /*! face area */
	  double & area()   { return m_area;   };

  };

  /*! \brief CViewerVertex class
  *   
  *   Vertex class for viewer
  *   Trait : vertex rgb color
  */
  class CCurvatureVertex : public  CVertex
  {
  protected:
	/*! curvature */
	double   m_curvature;
  public:
	  /*! CViewerVertex Constructor */
	  CCurvatureVertex()
	  {
		  m_curvature = 0; 
	  }
	  /*! vertex curvature */
	  double & k() { return m_curvature; };
  };
 

  /*! \brief CCurvatureEdge class
  *   
  *   Edge class for viewer
  *   Trait : Edge sharp
  */
  class CCurvatureEdge : public  CEdge
  {
  protected:
	/*! edge length */
	double   m_length;
	
  public:
	  /*! CViewerVertex Constructor */
	  CCurvatureEdge()
	  {
		  m_length = 0;
	  }

	 /*!
	  *	edge length
	  */
	 double& length() { return m_length; };

  };
 


  /*! \brief CCurvatureEdge class
  *   
  *   Edge class for viewer
  *   Trait : Edge sharp
  */
  class CCurvatureHalfEdge : public  CHalfEdge
  {
  protected:
	/*! corner angle */
	double   m_angle;
	
  public:
	  /*! CViewerVertex Constructor */
	  CCurvatureHalfEdge()
	  {
		  m_angle = 0;
	  }

	  /*!
	  *	corner angle
	  */
	 double& angle() { return m_angle; };

  };

/*-------------------------------------------------------------------------------------------------------------------------------------

	Viewer Mesh

--------------------------------------------------------------------------------------------------------------------------------------*/
/*! \brief CViewerMesh class
*
*	mesh class for viewer
*
*/
 template<typename V, typename E, typename F, typename H>
class CCurvatureMesh : public CBaseMesh<V,E,F,H>
{
public:
	
	typedef V CVertex;
	typedef E CEdge;
	typedef F CFace;
	typedef H CHalfEdge;

	typedef CBoundary<V,E,F,H> CBoundary;
	typedef CLoop<V,E,F,H> CLoop;
	typedef MeshVertexIterator<V,E,F,H> MeshVertexIterator;
	typedef MeshFaceIterator<V,E,F,H>   MeshFaceIterator;
	typedef MeshEdgeIterator<V,E,F,H>   MeshEdgeIterator;
	typedef VertexVertexIterator<V,E,F,H> VertexVertexIterator;
	typedef FaceVertexIterator<V,E,F,H> FaceVertexIterator;
	typedef VertexEdgeIterator<V,E,F,H> VertexEdgeIterator;
	typedef VertexFaceIterator<V,E,F,H> VertexFaceIterator;
	typedef VertexInHalfedgeIterator<V,E,F,H>  VertexInHalfedgeIterator;
	typedef FaceHalfedgeIterator<V,E,F,H> FaceHalfedgeIterator;
	
public:
};

typedef CCurvatureMesh<CCurvatureVertex, CCurvatureEdge, CCurvatureFace, CCurvatureHalfEdge> CKMesh;

unsigned long long CKMesh::m_input_traits  = 0;
unsigned long long CKMesh::m_output_traits = 0;

}
#endif  _CURVATURE_MESH_H_