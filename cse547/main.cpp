
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
void compute_edge_length( M * pMesh )
{
}

//compute corner angle

template<typename M>
void compute_corner_angle( M * pMesh )
{
}

//compute vertex curvature

template<typename M>
void compute_vertex_curvature( M * pMesh )
{
}

//verify Gauss-Bonnet theorem

template<typename M>
double total_curvature( M * pMesh )
{
	return 0;
}

//compute the Euler number

template<typename M>
int Euler_number( M * pMesh )
{
	return 0;
}


//compute the normals to each face

template<typename M>
void compute_face_normal( M* pMesh )
{
}

//compute the normals to each vertex

template<typename M>
void compute_vertex_normal( M* pMesh )
{
}

//compute the number of boundary loops
template<typename M>
int compute_boundary_loops( M * pMesh )
{
	return 0;
}


/*! main function for viewer
*/
int main( int argc, char * argv[] )
{
		CKMesh mesh;
		mesh.read_m( argv[1] );

		compute_face_normal<CKMesh>( &mesh );
		compute_vertex_normal<CKMesh>( &mesh );

		compute_edge_length<CKMesh>( &mesh );
		compute_corner_angle<CKMesh>( &mesh );
		compute_vertex_curvature<CKMesh>( &mesh );
		
		double K = total_curvature<CKMesh> ( &mesh );
		int E = Euler_number<CKMesh>( &mesh );
		int b = compute_boundary_loops( &mesh );

		int genus = (2-b-E)/2;

		std::cout << "Total curvature : " << K/(2*PI ) << " 2PI " << std::endl;
		std::cout << "Euler Number    : " << E << std::endl;
		std::cout << "Genus           : " << genus << std::endl;
		std::cout << "Boundary        : " << b << std::endl;

		std::cin.get();

		return 0;
}