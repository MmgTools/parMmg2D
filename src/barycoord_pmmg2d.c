/* =============================================================================
**  This file is part of the parmmg software package for parallel tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux, 2017-
**
**  parmmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  parmmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with parmmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the parmmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file barycoord_pmmg2d.c
 * \brief Barycentric coordinates for point localization in a mesh.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg2d.h"
#include "barycoord_pmmg2d.h"

/**
 * \param coord pointer to a double array
 * \param barycoord pointer to the point barycentric coordinates in the current triangle
 * \param ndim space dimension of the manifold
 *
 *  Get barycentric coordinates.
 *
 */
void PMMG2D_barycoord_get( double *val, PMMG2D_barycoord *phi, int ndim ) {
  int i;

  for( i = 0; i < ndim; i++ ) val[phi[i].idx] = phi[i].val;

}

int PMMG2D_barycoord_isInside( PMMG2D_barycoord *phi ) {
  if( phi[0].val > -PMMG2D_EPS )
    return 1;
  else
    return 0;
}

/**
 * \param mesh pointer to the mesh structure
 * \param ptr pointer to the current triangle
 * \param k index of the triangle
 * \param coord pointer to the point coordinates
 * \param normal unit normal of the current triangle
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if the point is outside the element, 1 if inside.
 *
 *  Compute the barycentric coordinates of a given point in a given triangle,
 *  sort them and evaluate if the point is inside.
 *
 */
int PMMG2D_barycoord2d_evaluate( MMG5_pMesh mesh,MMG5_pTria ptr,
                                 double *coord, PMMG2D_barycoord *barycoord ) {

  // Get barycentric coordinates and sort them in ascending order
  PMMG2D_barycoord2d_compute(mesh, ptr, coord, barycoord);
  qsort(barycoord,3,sizeof(PMMG2D_barycoord), PMMG2D_barycoord_compare);

  // Return inside/outside status 
  return PMMG2D_barycoord_isInside( barycoord );
}

/**
 * \param mesh pointer to the mesh structure
 * \param ptr pointer to the current triangle
 * \param k index of the triangle
 * \param coord pointer to the point coordinates
 * \param normal unit normal of the current triangle
 * \param barycoord pointer to the point barycentric coordinates in the current
 * tetra
 *
 * \return 0 if fail, 1 if success
 *
 *  Compute the barycentric coordinates of a given point in a given triangle.
 *
 */
int PMMG2D_barycoord2d_compute( MMG5_pMesh mesh,MMG5_pTria ptr,double *coord,
                                PMMG2D_barycoord *barycoord ) {
  double *c1,*c2,vol;
  int    ia;

  // Retrieve tria area
  vol = PMMG2D_quickarea(mesh->point[ptr->v[0]].c, mesh->point[ptr->v[1]].c, mesh->point[ptr->v[2]].c);

  // Retrieve face areas and compute barycentric coordinates
  for( ia = 0; ia < 3; ia++ ) {
    c1 = mesh->point[ptr->v[MMG5_inxt2[ia]]].c;
    c2 = mesh->point[ptr->v[MMG5_inxt2[ia+1]]].c;
    barycoord[ia].val = PMMG2D_quickarea( coord, c1, c2 )/vol;
    barycoord[ia].idx = ia;
  }

  return 1;
}

/**
 * \param a pointer to point barycentric coordinates
 * \param b pointer to point barycentric coordinates
 *
 * \return -1 if (a < b), +1 if (a > b), 0 if equal
 *
 *  Compare the barycentric coordinates of a given point in a given triangle.
 *
 */
int PMMG2D_barycoord_compare( const void *a,const void *b ) {
  PMMG2D_barycoord *coord_a;
  PMMG2D_barycoord *coord_b;

  coord_a = (PMMG2D_barycoord *)a;
  coord_b = (PMMG2D_barycoord *)b;

  if( coord_a->val > coord_b-> val ) return 1;
  if( coord_a->val < coord_b-> val ) return -1;

  return 0;
}

int PMMG2D_barycoord_isBorder( PMMG2D_barycoord *phi,int *ifoundEdge,int *ifoundVertex ) {
  if( phi[0].val < PMMG2D_EPS ) {
    if( phi[1].val < PMMG2D_EPS ) {
      *ifoundVertex = phi[2].idx;
    }
    else {
      *ifoundEdge = phi[0].idx;
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param k index of the triangle to analyze
 * \param ppt pointer to the point to locate
 * \param barycoord barycentric coordinates of the point to be located
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a surface point in its closest background triangles, and find its
 *  closest point.
 *
 */
int PMMG2D_barycoord2d_getClosest( MMG5_pMesh mesh,int k,MMG5_pPoint ppt,
                                   PMMG2D_barycoord *barycoord ) {
  MMG5_pTria ptr;
  double *c,dist[2],norm,min;
  int i,d,itarget;

  ptr = &mesh->tria[k];

  c = mesh->point[ptr->v[0]].c;
  for( d = 0; d < 2; d++ )
    dist[d] = ppt->c[d] - c[d];
  norm = sqrt(dist[0]*dist[0]+dist[1]*dist[1]);
  min = norm;
  itarget = 0;

  for( i = 1; i < 3; i++ ) {
    c = mesh->point[ptr->v[i]].c;
    for( d = 0; d < 2; d++ )
      dist[d] = ppt->c[d] - c[d];
    norm = sqrt(dist[0]*dist[0]+dist[1]*dist[1]);
    if( norm < min ) {
      min = norm;
      itarget = i;
    }
  }

  for( i = 0; i < 3; i++ ) {
    barycoord[i].val = 0.0;
    barycoord[i].idx = i;
  }
  barycoord[itarget].val = 1.0;

  return 1;
}

/**
 * \param a coordinates of the first point
 * \param b coordinates of the second point
 * \param c coordinates of the third point
 *
 * \return the area of the triangle
 *
 */
double PMMG2D_quickarea(double *a,double *b,double *c)
{return 0.5*((b[0] - a[0])*(c[1] - a[1]) - (b[1] - a[1])*(c[0] - a[0]));}
