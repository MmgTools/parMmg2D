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
 * \file locate_pmmg2d.h
 * \brief Point localization for interpolation on a new mesh.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg2d.h"
#include "locate_pmmg2d.h"

/**
 * \param mesh pointer to the background mesh structure
 * \param ptr pointer to the triangle to analyze
 * \param k index of the triangle
 * \param ppt pointer to the point to locate
 * \param triaNormal unit normal of the current triangle
 * \param barycoord barycentric coordinates of the point to be located
 * \param h distance from the triangle
 * \param closestDist distance from the closest triangle
 * \param closestTria index of the closest triangle (with negative sign)
 *
 * \return 1 if found; 0 if not found
 *
 *  Locate a point in a background triangle, and provide its barycentric
 *  coordinates.
 *
 */
int PMMG2D_locatePointInTria( MMG5_pMesh mesh,MMG5_pTria ptr,int k,MMG5_pPoint ppt,
                              PMMG2D_barycoord *barycoord,
                              double *closestDist,int *closestTria ) {
  MMG5_pPoint    ppt0;
  double         norm,dist[2];
  int            j,d,found;

  // Mark tria 
  ptr->flag = mesh->base;

  // Evaluate point in triangle through barycentric coordinates 
  found = PMMG2D_barycoord2d_evaluate( mesh,ptr,ppt->c,barycoord );

  // Distance from center of mass 
  for( d = 0; d < 2; d++ )
    dist[d] = ppt->c[d];
  for( j = 0; j < 3; j++ ) {
    ppt0 = &mesh->point[ptr->v[j]];
    for( d = 0; d < 2; d++ )
      dist[d] -= ppt0->c[d]/3.0;
  }
  norm = 0;
  for( d = 0; d < 2; d++ )
    norm += dist[d]*dist[d];
  norm = sqrt(norm);

  // Save element index if it is the closest one 
  if( norm < *closestDist ) {
    *closestDist = norm;
    *closestTria = k;
  }
  assert(*closestTria);

  return found;
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param triaNormals non-normalized triangle normals of all mesh triangles
 * \param iTria pointer to the index of the found triangle
 * \param closestTria pointer to the index of the closest triangle
 * \param closestDist pointer to the closest distance
 *
 * \return 1 if found, 0 otherwise.
 *
 *  Exhaustive point search on the background triangles. If the point is not
 *  found, the triangle pointers points to the closest triangle.
 *
 */
int PMMG2D_locatePoint_exhaustTria( MMG5_pMesh mesh,MMG5_pPoint ppt,
                                    PMMG2D_barycoord *barycoord,
                                    int *iTria,int *closestTria,double *closestDist ) {
  MMG5_pTria     ptr;

  for( *iTria = 1; *iTria <= mesh->nt; (*iTria)++ ) {

    // Increase step counter
    ppt->s--;

    // Get tetra
    ptr = &mesh->tria[*iTria];
    if ( !MG_EOK(ptr) ) continue;

    //¨Skip already analized tetras
    if( ptr->flag == mesh->base ) continue;

    // Exit the loop if you find the element
    if( PMMG2D_locatePointInTria( mesh, ptr, *iTria, ppt,
                                  barycoord,
                                  closestDist, closestTria ) ) break;

  }

  if( *iTria <= mesh->nt ) {
    return 1;
  }
  else {
    *iTria = *closestTria;
    // Recompute barycentric coordinates
    if( !PMMG2D_locatePointInTria( mesh, ptr, *iTria, ppt,
                                   barycoord,
                                   closestDist, closestTria ) ) {
      // Recompute barycentric coordinates to the closest point
      PMMG2D_barycoord2d_getClosest( mesh,*iTria,ppt,barycoord );
    }
    return 0;
  }
}

/**
 * \param mesh pointer to the background mesh structure
 * \param ppt pointer to the point to locate
 * \param triaNormals unit normals of the all triangles in the mesh
 * \param barycoord barycentric coordinates of the point to be located
 * \param iTria pointer to the index of the triangle
 * \param ifoundEdge pointer to the index of the local edge
 * \param ifoundVertex pointer to the index of the local vertex
 *
 * \return 0 if not found (closest), 1 if found, -1 if found through exhaustive
 * search.
 *
 *  Locate a point in a background mesh surface by traveling the triangles
 *  adjacency.
 *
 */
int PMMG2D_locatePoint( MMG5_pMesh mesh, int* list_triangles, MMG5_pPoint ppt,
                        PMMG2D_barycoord *barycoord,
                        int *iTria,int *ifoundEdge,int *ifoundVertex ) {
  MMG5_pTria     ptr;
  int            i, k, closestTria;
  double         closestDist;
  static int     mmgWarn0=0, mmgWarn1=0;
  int            ier;
  int            take_closest = 1;

  ++mesh->base;

  closestTria = 0;
  closestDist = 1.0e10;

  *ifoundEdge   = PMMG2D_UNSET;
  *ifoundVertex = PMMG2D_UNSET;

  for (i = 1; i <= list_triangles[0]; i++) {

    k = list_triangles[i];
    ptr = &mesh->tria[k];
    if ( !MG_EOK(ptr) ) continue;

    // Exit the loop if you find the element
    if( PMMG2D_locatePointInTria( mesh, ptr, k, ppt,
                                  barycoord, &closestDist, &closestTria ) ) {
      PMMG2D_barycoord_isBorder( barycoord, ifoundEdge, ifoundVertex );
      *iTria = k;
      return 1;
    }

  }

  // Recompute barycentric coordinates to the closest point
  *iTria = closestTria;
  PMMG2D_barycoord2d_getClosest( mesh,*iTria,ppt,barycoord );
  if (take_closest) return 0;

  if ( !mmgWarn0 ) {
    mmgWarn0 = 1;
    if ( mesh->info.imprim > PMMG2D_VERB_DETQUAL ) {
      fprintf(stderr,"\n  ## Warning %s: Cannot locate point,"
              " performing exhaustive research.\n",__func__);
    }
  }

  ier = PMMG2D_locatePoint_exhaustTria( mesh, ppt, barycoord,
                                        iTria, &closestTria, &closestDist );

  if( ier ) {
    return -1;
  } 
  else {
    // Element not found: Return the closest one 
    if ( !mmgWarn1 ) {
      mmgWarn1 = 1;
      if ( mesh->info.imprim > PMMG2D_VERB_VERSION ) {
        fprintf(stderr,"\n  ## Warning %s: Point not located, smallest external area %e.",
                __func__,closestDist);
      }
    }
    return 0;
  }

}

/**
 * \param mesh pointer to the current mesh structure
 * \param ip point index
 * \param ier error code
 * \param myrank process rank
 * \param igrp mesh group index
 *
 * Analise the found element and display warnings for exhaustive searches and
 * points not found.
 *
 */
void PMMG2D_locatePoint_errorCheck( MMG5_pMesh mesh, int ip, int ier,
                                   int myrank) {
  MMG5_pPoint ppt;
  static int8_t pmmgWarn0 = 0;
  static int8_t pmmgWarn1 = 0;

  ppt = &mesh->point[ip];

  if( !ier ) {
    if ( !pmmgWarn0 ) {
      pmmgWarn0 = 1;
      fprintf(stderr,"\n  ## Warning: %s (rank %d): at least one"
              " localisation issue: closest element for"
              " point %d (tag %d), coords %e %e %e\n",__func__,myrank,
              ip,ppt->tag,ppt->c[0],ppt->c[1],ppt->c[2]);
    }
  } else if ( ier < 0 ) {
    if ( !pmmgWarn1 ) {
      pmmgWarn1 = 1;
      fprintf(stderr,"\n  ## Warning: %s (rank %d): at least one"
              " exhaustive search for"
              " point %d (tag %d), coords %e %e %e\n",__func__,myrank,
              ip,ppt->tag,ppt->c[0],ppt->c[1],ppt->c[2]);
    }
  }
}
