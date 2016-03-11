
/*******************************************************************************
*
*      CSE 564 Visualization @ Stony Brook University
*
*      David Gu Feb 16, 2015
*
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CurvatureMesh/CurvatureMesh.h"



using namespace MeshLib;


//compute edge length

template<typename M>
void compute_edge_length(M * pMesh)
{
	for (M::MeshEdgeIterator eiter(pMesh); !eiter.end(); eiter++)
	{
		M::CEdge * pEdge = *eiter;
		//M::CHalfEdge * ph = pEdge;
		CVertex * s_v = pEdge->halfedge(0)->source();
		CVertex * t_v = pEdge->halfedge(0)->target();
		pEdge->length() = (t_v->point() - s_v->point()).norm();
	}
}

//compute corner angle

template<typename M>
void compute_corner_angle(M * pMesh)
{
	for (M::MeshFaceIterator fiter(pMesh); !fiter.end(); fiter++)
	{
		M::CFace * pF = *fiter;
		for (M::FaceHalfedgeIterator hiter(pF); !hiter.end(); ++hiter)
		{
			M::CHalfEdge * pH = *hiter;
			/*M::CEdge * pEdge = pH->edge();
			double li = pH->length();
			double lj = pH->he_prev()->edge()->length();
			double lk = pH->he_next()->edge()->length();*/
			CVertex * v0 = pH->source();
			CVertex * v1 = pH->he_next()->source();
			CVertex * v2 = pH->he_prev()->source();
			double li = (v1->point() - v0->point()).norm();
			double lj = (v2->point() - v1->point()).norm();
			double lk = (v0->point() - v2->point()).norm();
			double theta = acos((lj*lj + lk*lk - li*li) / (2 * lj*lk));
			//store the corner angle at pH->angle();
			pH->angle() = theta;
		}
	}
}

//compute vertex curvature

template<typename M>
void compute_vertex_curvature(M * pMesh)
{
	for (M::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
	{
		M::CVertex * pV = *viter;
		//Fill in your Code
		double cur_sum = 0.0;
		for (M::VertexInHalfedgeIterator vhiter(pMesh, pV); !vhiter.end(); ++vhiter)
		{
			M::CHalfEdge * pH = *vhiter;
			cur_sum += pH->angle();
		}
		// curvature stores in pV->k();
		if (pMesh->isBoundary(pV)){
			pV->k() = PI - cur_sum;
		}
		else{
			pV->k() = 2 * PI - cur_sum;
		}

	}
}

//verify Gauss-Bonnet theorem
/*int from_Euler =2* PI * Euler_number(pMesh);
	int from_Gauss = total_curvature(pMesh);

	return (from_Euler-from_Gauss);*/

template<typename M>
double total_curvature(M * pMesh)
{
	double result = 0.0;
	for (M::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
	{
		M::CVertex * pV = *viter;
		result += pV->k();
	}
	return result;
}

//compute the Euler number

template<typename M>
int Euler_number(M * pMesh)
{
	int f_count = 0;
	int e_count = 0;
	int v_count = 0;
	for (M::MeshFaceIterator fiter(pMesh); !fiter.end(); fiter++)
	{
		f_count++;
	}
	for (M::MeshEdgeIterator eiter(pMesh); !eiter.end(); eiter++)
	{
		e_count++;
	}
	for (M::MeshVertexIterator viter(pMesh); !viter.end(); viter++)
	{
		v_count++;
	}

	int result = v_count + f_count - e_count;
	return result;
}


//compute the normals to each face

template<typename M>
void compute_face_normal(M* pMesh)
{
	for (M::MeshFaceIterator fiter(pMesh); !fiter.end(); fiter++)
	{
		M::CFace * pFace = *fiter;
		CVertex * v_arr[3];
		int i = 0;
		for (M::FaceVertexIterator fviter(pFace); !fviter.end(); ++fviter)
		{
			CVertex * pV = *fviter;
			v_arr[i++] = pV;
			
		}
		CPoint d = (v_arr[1]->point() - v_arr[0]->point()) ^ (v_arr[2]->point() - v_arr[0]->point());
		//area is stored in pFace->area();
		pFace->area() = d.norm();
		//normal is stored in pFace->normal();
		pFace->normal() = d / d.norm();
	}
}

//compute the normals to each vertex

template<typename M>
void compute_vertex_normal(M* pMesh)
{
	for (M::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
	{
		M::CVertex * pV = *viter;
		int i = 0;
		double farea_sum = 0;
		CPoint fnormal_sum(0.0, 0.0, 0.0);
		for (M::VertexFaceIterator vfiter(pV); !vfiter.end(); ++vfiter)
		{
			M::CFace * pF = *vfiter;
			fnormal_sum += pF->normal() * pF->area();
			farea_sum += pF->area();
			i++;
		}
		assert(farea_sum > 0.0);
		CPoint d = fnormal_sum / farea_sum;

		//vertex normal is stored in pV->normal();
		pV->normal() = d / d.norm();
	}

}

//compute the number of boundary loops
template<typename M>
int compute_boundary_loops(M * pMesh)
{
	int result = 0;
	
	std::list<CHalfEdge*> boundary;
	std::vector<int> boundary_length;
	int * edge_flag = new int[pMesh->numEdges()];
	int cc = 0; // for debug
	for (M::MeshEdgeIterator eiter(pMesh); !eiter.end(); eiter++)
	{
		cc++;
		int bl = 0; // boundary length
		CEdge * edge = *eiter;
		//he = M::edgeHalfedge(edge, 1);
		if ((edge->halfedge(1) == NULL)&&(*(edge_flag + edge->id()-1) != 1))
		{
			
			boundary.push_back(edge->halfedge(0));
			bl++;
			*(edge_flag + edge->id() - 1) = 1;
			CHalfEdge * he = edge->halfedge(0);
			while (true)
			{

				if (bl>1 && he == edge->halfedge(0))
				{
					boundary_length.push_back(boundary.size());
					break;
				}
				else if (he->he_next()->he_sym() == NULL)
				{
					boundary.push_back(he->he_next());
					bl++;
					*(edge_flag + he->he_next()->edge()->id() - 1) = 1;
					he = he->he_next();
				}
				else
				{
					he = he->he_next()->he_sym();
				}
			}
		}
	}
	delete [] edge_flag;
	result = boundary_length.size();
	return result;
} 


/*! main function for viewer
*/
int main(int argc, char * argv[])
{
	CKMesh mesh;
	mesh.read_m(argv[1]);

	compute_face_normal<CKMesh>(&mesh);
	compute_vertex_normal<CKMesh>(&mesh);

	compute_edge_length<CKMesh>(&mesh);
	compute_corner_angle<CKMesh>(&mesh);
	compute_vertex_curvature<CKMesh>(&mesh);

	double K = total_curvature<CKMesh>(&mesh);
	int E = Euler_number<CKMesh>(&mesh);
	int b = compute_boundary_loops(&mesh);

	int genus = (2 - b - E) / 2;

	std::cout << "Total curvature : " << K / (2 * PI) << " 2PI " << std::endl;
	std::cout << "Euler Number    : " << E << std::endl;
	std::cout << "Genus           : " << genus << std::endl;
	std::cout << "Boundary        : " << b << std::endl;

	std::cin.get();

	return 0;
}
