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
 * \file interpmesh_pmmg2d.c
 * \brief Interpolate data from a background mesh to the current mesh.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "parmmg2d.h"
#include "interpmesh_pmmg2d.h"

/**
 * \param m initial symmetric matrix.
 * \param mi inverted matrix.
 *
 * Invert 2x2 non-symmetric matrix stored in 2 dimensions
 *
 */
static int PMMG2D_invmat22(double *m,double *mi) {
  double det = m[0]*m[2] - m[1]*m[1];
  if ( fabs(det) < 1.e-16 ) return 0;

  det = 1/det;
  mi[0] = det*m[2];
  mi[1] = -det*m[1];
  mi[2] = det*m[0];

  return 1;
}


/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param l local index of the edge on the background triangle
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background edge.
 */
int PMMG2D_interp2bar_iso( MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol oldMet,
                           MMG5_pTria ptr, int ip, int l, PMMG2D_barycoord *barycoord ) {
  int    i0,i1;
  double phi[3];

  assert( met->size == 1 );

  PMMG2D_barycoord_get( phi, barycoord, 3 );

  i0 = MMG5_inxt2[l];
  i1 = MMG5_iprv2[l];

  // Linear interpolation of the squared size
  met->m[ip] = phi[i0]*oldMet->m[ptr->v[i0]] +
               phi[i1]*oldMet->m[ptr->v[i1]];

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param l local index of the edge on the background triangle
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point the metrics inverse on a target background
 *  edge.
 *
 */
int PMMG2D_interp2bar_ani( MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol oldMet,
                           MMG5_pTria ptr, int ip, int l, PMMG2D_barycoord *barycoord ) {
  int    i0,i1;
  double phi[3],mi[2][3],mint[3];
  int    isize,nsize;

  assert( met->size == 3 );
  nsize  = met->size;

  PMMG2D_barycoord_get( phi, barycoord, 3 );

  i0 = MMG5_inxt2[l];
  i1 = MMG5_iprv2[l];

  if( !PMMG2D_invmat22( &oldMet->m[nsize*ptr->v[i0]], mi[0] ) ) return 0;
  if( !PMMG2D_invmat22( &oldMet->m[nsize*ptr->v[i1]], mi[1] ) ) return 0;

  // Linear interpolation of the metrics
  for( isize = 0; isize < nsize; isize++ ) {
    mint[isize] = phi[i0]*mi[0][isize]+
                  phi[i1]*mi[1][isize];
  }

  if( !PMMG2D_invmat22( mint, &met->m[nsize*ip] ) ) return 0;

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point metrics on a target background triangle.
 *  This function is analogous to the MMG5_interp2bar_iso() function.
 */
int PMMG2D_interp3bar_iso( MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol oldMet,
                           MMG5_pTria ptr, int ip, PMMG2D_barycoord *barycoord ) {
  double phi[3];
  int iadr,i,j;

  iadr = met->size *ip;

  assert (mesh->npmax == met->npmax );

  PMMG2D_barycoord_get( phi, barycoord, 3 );

  // Linear interpolation of the squared size
  for ( j=0; j<met->size; ++j ) {
    met->m [ iadr + j ] = 0.0;
  }

  for( i=0; i<3; i++ ) {
    for ( j=0; j<met->size; ++j ) {
      // Barycentric coordinates could be permuted
      met->m[iadr+j] += phi[i]*oldMet->m[ptr->v[i]*met->size+j];
    }
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param ptr pointer to the target background triangle
 * \param ip index of the current point
 * \param barycoord barycentric coordinates of the point to be interpolated
 *
 * \return 0 if fail, 1 if success
 *
 *  Linearly interpolate point the metrics inverse on a target background
 *  triangle.
 *  This function is analogous to the MMG5_interp4barintern() function.
 *
 */
int PMMG2D_interp3bar_ani( MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol oldMet,
                           MMG5_pTria ptr, int ip, PMMG2D_barycoord *barycoord ) {
  double phi[3],mi[3][3],mint[3];
  int    i,isize,nsize;

  assert( met->size == 3 );
  nsize  = met->size;

  PMMG2D_barycoord_get( phi, barycoord, 3 );

  for( i=0; i<3; i++ ) {
    if( !PMMG2D_invmat22( &oldMet->m[nsize*ptr->v[i]], mi[i]) ) return 0;
  }

  // Linear interpolation of the metrics
  for( isize = 0; isize < nsize; isize++ ) {
    mint[isize] = phi[0]*mi[0][isize]+
                  phi[1]*mi[1][isize]+
                  phi[2]*mi[2][isize];
  }

  if( !PMMG2D_invmat22( mint, &met->m[nsize*ip] ) ) return 0;

  return 1;
}

/**
 * \param oldMesh pointer to the background mesh
 * \param met pointer to the current metrics
 * \param oldMet pointer to the background metrics
 * \param c1 x-coordinate of the point to be interpolated
 * \param c2 y-coordinate of the point to be interpolated
 * \param ifoundTria1 index to the background triangle containing the point
 * \param ip index of the current point
 *
 * \return 0 if fail, 1 if success
 *
 *  Barycentric interpolation of the metrics on a target background
 *  list of points, based on their inverse distance to the point.
 *
 */
int PMMG2D_extended_interpNbar( MMG5_pMesh oldMesh, MMG5_pSol met, MMG5_pSol oldMet, 
                                int nlayers, double c1, double c2, int ifoundTria1, int ip )
{
  int i, j, k, l, m, n, p;
  int k1, nsize;
  double distance, sum_distance;
  int *list_triangles = (int*)malloc(sizeof(int));
  int *list_points = (int*)malloc(0*sizeof(int));

  nsize = met->size;
  for (k = 1; k <= oldMesh->nt; k++) oldMesh->tria[k].flag = 0;
  oldMesh->base++;

  // 1) Compute list of adjacent triangles which will participate to the interpolation
  list_triangles[0] = ifoundTria1;
  l = 1;
  for (i = 0; i < nlayers; i++) {

    int size = l;
    int list_initial_triangles[size];
    for (n = 0; n < size; n++) list_initial_triangles[n] = list_triangles[n];

    for (j = 0; j < size; j++) { // Loop over the list of triangles at the step i

      k = list_initial_triangles[j];

      for( n = 1; n < 4; n++ ) { // Loop over the adjacent triangles to triangle k

        k1 = oldMesh->adja[3*(k-1)+n];
        if( !k1 ) continue;
        k1 /= 3;
  
        if ( (&oldMesh->tria[k1])->flag ) continue; // The triangle is already exchanged
  
        oldMesh->tria[k1].flag = oldMesh->base;
        
        l++; // Number of triangles in the list
  
        list_triangles = (int*)realloc(list_triangles, l*sizeof(int));
        list_triangles[l-1] = k1;
      }

    }

  }

  // 2) Store the points participating to the interpolation once
  for (k = 1; k <= oldMesh->np; k++) oldMesh->point[k].flag = 0;

  p = 0;
  for (i = 0; i < l; i++) {
    k = list_triangles[i];
    for (j = 0; j < 3; j++) {
      n = oldMesh->tria[k].v[j];

      if ( (&oldMesh->point[n])->flag ) continue;

      oldMesh->point[n].flag = oldMesh->base;

      p++; // Number of points in the list
      list_points = (int*)realloc(list_points, p*sizeof(int));
      list_points[p-1] = n;
    }
  }

  free(list_triangles);

  // 3) Compute the inverse distance to each point and perform the interpolation
  for (i = 0; i < nsize; i++) met->m[nsize*ip+i] = 0.;
  sum_distance = 0.;

  for (i = 0; i < p; i++) {
    n = list_points[i];
    distance = sqrt( (c1-oldMesh->point[n].c[0]) * (c1-oldMesh->point[n].c[0]) +
                     (c2-oldMesh->point[n].c[1]) * (c2-oldMesh->point[n].c[1]) );

    if (distance < PMMG2D_EPS) {
      for (j = 0; j < nsize; j++) met->m[nsize*ip+j] = oldMet->m[nsize*n+j];
      return 1;
    }

    sum_distance += 1. / distance;
    for (j = 0; j < nsize; j++) {
      met->m[nsize*ip+j] += oldMet->m[nsize*n+j] / distance;
    }
  }

  for (j = 0; j < nsize; j++) met->m[nsize*ip+j] /= sum_distance;

  free(list_points);

  return 1;
}

/**
 * \param mesh pointer to the current mesh.
 * \param met pointer to the current metrics.
 * \param oldMesh pointer to the background mesh.
 * \param oldMet pointer to the background metrics.
 * \param idest index of the target point.
 * \param isrc index of the source point.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric of a point.
 *
 */
int PMMG2D_copyMetrics( MMG5_pMesh mesh, MMG5_pSol met, MMG5_pMesh oldMesh,
                        MMG5_pSol oldMet, int idest, int isrc ) {
  int isize,nsize;

  nsize = met->size;

  for( isize = 0; isize<nsize; isize++ ) {
    met->m[nsize*idest+isize] = oldMet->m[nsize*isrc+isize];
  }

  return 1;
}

/**
 * \param mesh pointer to the current mesh.
 * \param oldMesh pointer to the background mesh.
 * \param met pointer to the current metrics.
 * \param oldMet pointer to the background metrics.
 * \param permNodGlob permutation array for nodes.
 * \param inputMet 1 if user provided metric.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the data stored in a solution structure of a freezed interface point.
 *
 */
static
int PMMG2D_copySol_point( MMG5_pMesh mesh, MMG5_pMesh oldMesh,
                          MMG5_pSol sol, MMG5_pSol oldSol, int npt ) {
  MMG5_pPoint    ppt;
  int nsize, i, ip;
  int* metric_req = (int*) malloc(npt * sizeof(int));

  nsize   = sol->size;

  // 1) Store the old metric values of the required edges
  for(ip = 1; ip <= oldMesh->np; ip++) {
    ppt = &oldMesh->point[ip];
    if( !MG_VOK(ppt) ) continue;
    if( (ppt->tag & MG_REQ) || (ppt->tag & MG_PARBDY) ) metric_req[ppt->ref] = ip;
  }

  // 2) Copy the metric values
  for (ip = 1; ip <= mesh->np; ip++) {
    ppt = &mesh->point[ip];
    if( !MG_VOK(ppt) ) continue;
    if( (ppt->tag & MG_REQ) || (ppt->tag & MG_PARBDY) ) {
      for( i = 0; i < nsize; i++ ) {
        sol->m[nsize*ip+i] = oldSol->m[nsize*metric_req[ppt->ref] + i];
      }
    }
  }

  free(metric_req);

  return 1;
}

/**
 * \param mesh pointer to the current mesh.
 * \param oldMesh pointer to the background mesh.
 * \param met pointer to the current metrics.
 * \param oldMet pointer to the background metrics.
 * \param permNodGlob permutation array for nodes.
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the metric of a freezed interface point.
 *
 */
int PMMG2D_copyMetrics_point( MMG5_pMesh mesh, MMG5_pMesh oldMesh,
                              MMG5_pSol met, MMG5_pSol oldMet, int npt ) {
  int ier;

  if ( mesh->info.hsiz > 0.0 ) return 1;

  ier =  PMMG2D_copySol_point( mesh, oldMesh, met, oldMet, npt );

  return ier;
}

static int PMMG2D_is_point_belong_to_line( double x1, double y1, double x2, double y2, double px, double py ) {
  int r = fabs((x2 - x1)*(py - y1) - (px - x1)*(y2 - y1)) < PMMG2D_EPS;
  if (!r) return 0;

  if (fabs(x1 - x2) > PMMG2D_EPS) {
    r = ( x1 <= px && px <= x2 ) || ( x1 >= px && px >= x2 );
  }
  else {
    r = ( y1 <= py && py <= y2 ) || ( y1 >= py && py >= y2 );
  }
  return r;
}

/**
 * \param size_polygon is the number of points in the polygon.
 * \param polygon contains the points and their coordinates.
 * \param px is the x-coordinate of the point.
 * \param py is the y-coordinate of the point.
 *
 * \return the winding number
 *
 * Compute the winding number associated to the points (px,py) and the given polygon
 *
 */
static int PMMG2D_winding_number(int size_polygon, double** polygon, double px, double py) {

  int winding_number = 0;
  double x1, x2, y1, y2;
  int i;

  for (i = 0; i < size_polygon; i++) {
    x1 = polygon[i][0];
    y1 = polygon[i][1];
    x2 = polygon[(i+1)%size_polygon][0];
    y2 = polygon[(i+1)%size_polygon][1];

    if (PMMG2D_is_point_belong_to_line(x1,y1,x2,y2,px,py)) return 1;

    if (y1 <= py) {
      if (y2 <= py) continue;
      if ((x2 - x1)*(py - y1) <= (px - x1)*(y2 - y1)) continue;
      winding_number++;
    }
    else {
      if (y2 > py) continue;
      if ((x2 - x1)*(py - y1) > (px - x1)*(y2 - y1)) continue;
      winding_number--;
    }
  }

  return winding_number;

}

/**
 * \param size_polygon is the number of points in the polygon.
 * \param polygon contains the points and their coordinates.
 * \param px is the x-coordinate of the point.
 * \param py is the y-coordinate of the point.
 *
 * \return 1 if the point is inside the polygon, 0 otherwise
 *
 */
static int PMMG2D_is_point_inside_polygon(int size_polygon, double** polygon, double px, double py) { 
  if (PMMG2D_winding_number(size_polygon, polygon, px, py)) {
    return 1;
  }
  else return 0;
}

/**
 * \param mesh pointer to the current mesh structure.
 * 
 * Compute the coordinates (xmin, xmax, ymin, ymax) of the smallest box containing every point in the mesh
 *
*/
static void find_bounding_box(MMG5_pMesh mesh, double* xmin, double* xmax, double* ymin, double* ymax) {
  *xmax = -1.e20;
  *xmin = 1.e20;
  *ymax = -1.e20;
  *ymin = 1.e20;

  for (int k = 1; k <= mesh->np; k++) {

    if (mesh->point[k].c[0] < *xmin) *xmin = mesh->point[k].c[0];
    if (mesh->point[k].c[0] > *xmax) *xmax = mesh->point[k].c[0];
    if (mesh->point[k].c[1] < *ymin) *ymin = mesh->point[k].c[1];
    if (mesh->point[k].c[1] > *ymax) *ymax = mesh->point[k].c[1];

  }
}

/**
 * \param mesh pointer to the current mesh structure.
 * \param oldMesh pointer to the background mesh structure.
 * \param met pointer to the current metrics structure.
 * \param oldMet pointer to the background metrics structure.
 * \param faceAreas pointer to the array of oriented face areas.
 * \param triaNormals pointer to the array of non-normalized triangle normals.
 * \param permNodGlob permutation array of nodes.
 * \param inputMet 1 if user provided metric.
 * \param myrank process rank.
 * \param igrp current mesh group.
 * \param locStats pointer to the localization statistics structure.
 *
 * \return 0 if fail, 1 if success
 *
 * Interpolate metrics and solution fields for all groups from background
 * to current meshes.
 * For the metric:
 *  Do nothing if no metrics is provided (inputMet == 0), otherwise:
 *  - if the metrics is constant, recompute it;
 *  - else, interpolate the non-constant metrics.
 *
 * For the solution fields: Do nothing if no solution field is provided
 *   (mesh->nsols == 0), interpolate the non-constant field otherwise.
 *
 *
 */
static
int PMMG2D_interpMetrics_mesh( MMG5_pMesh mesh, MMG5_pMesh oldMesh, MMG5_pMesh initialMesh,
                               MMG5_pSol met,MMG5_pSol oldMet, MMG5_pSol initialMet,
                               double ** polygon, int size_polygon,
                               uint8_t inputMet, int optim_interp, int nlayers, int myrank ) {
  MMG5_pTria pt;
  MMG5_pPoint ppt;
  PMMG2D_barycoord barycoord[3];
  int         ifoundTria1, ifoundTria2;
  int         ifoundEdge1, ifoundEdge2;
  int         ifoundVertex1, ifoundVertex2;
  int         ip, ie, iloc, idx, idy, i, k;
  int         ismet,ier;
  int***      list_triangles;
  int***      list_triangles2;
  int         GRID_SIZE = sqrt(mesh->nt); 
  double      xmin, xmax, ymin, ymax; 
  double      xmin2, xmax2, ymin2, ymax2;

  ismet = 1;
  if( inputMet != 1 ) {
    ismet = 0;
  }
  else  if( mesh->info.hsiz > 0.0 ) {
    // Compute constant metrics 
    if ( !MMG2D_Set_constantSize(mesh,met) ) return 0;
    ismet = 0;
  }

  if ( !ismet ) {
    // Nothing to do (no metric)
    return 1;
  }

  // Interpolate metrics
  if (optim_interp) {
    initialMesh->base = 0;
    for ( ip = 1; ip <= initialMesh->nt; ip++ ) {
      pt = &initialMesh->tria[ip];
      if ( !MG_EOK(pt) ) continue;
      pt->flag = initialMesh->base;
    }
    find_bounding_box(initialMesh, &xmin2, &xmax2, &ymin2, &ymax2);
    list_triangles2 = grid_size_triangles(initialMesh, xmin2, ymin2, xmax2, ymax2, GRID_SIZE);
    find_bounding_box(oldMesh, &xmin, &xmax, &ymin, &ymax);
    list_triangles = grid_size_triangles(oldMesh, xmin, ymin, xmax, ymax, GRID_SIZE);
  }
  else {
    oldMesh->base = 0;
    for ( ip = 1; ip <= oldMesh->nt; ip++ ) {
      pt = &oldMesh->tria[ip];
      if ( !MG_EOK(pt) ) continue;
      pt->flag = oldMesh->base;
    }
    find_bounding_box(oldMesh, &xmin, &xmax, &ymin, &ymax);
    list_triangles = grid_size_triangles(oldMesh, xmin, ymin, xmax, ymax, GRID_SIZE);
  }

  ifoundTria1 = 1;
  ifoundTria2 = 1;

  for( ie = 1; ie <= mesh->np; ie++ ) {
    ppt = &mesh->point[ie];
    ppt->flag = mesh->base;
  }

  // Loop on new triangles, and localize their vertices in the old mesh
  mesh->base++;
  for( ie = 1; ie <= mesh->nt; ie++ ) {
    pt = &mesh->tria[ie];

    if( !MG_EOK(pt) ) continue;

    for( iloc = 0; iloc < 3; iloc++ ) {
      ip = pt->v[iloc];
      ppt = &mesh->point[ip];

      if( !MG_VOK(ppt) ) continue;

      /* Skip already interpolated points */
      if( ppt->flag == mesh->base ) continue;

      if( (ppt->tag & MG_REQ) || (ppt->tag & MG_PARBDY) ) {
        /* Flag point as interpolated */
        ppt->flag = mesh->base;
        continue; // treated by copyMetric_points
      } 
      else {

        if (optim_interp && PMMG2D_is_point_inside_polygon(size_polygon, polygon, ppt->c[0], ppt->c[1]) ) {

          // Find the cell containing the point
          idx = (int)((ppt->c[0] - xmin2) / ((xmax2 - xmin2)/GRID_SIZE));
          idy = (int)((ppt->c[1] - ymin2) / ((ymax2 - ymin2)/GRID_SIZE));
          if (idx == -1) idx = 0; // Possible if the point is on the domaine boundary
          if (idx == GRID_SIZE) idx = GRID_SIZE-1;
          if (idy == -1) idy = 0;
          if (idy == GRID_SIZE) idy = GRID_SIZE-1;

          // Locate point in the initial mesh 
          ier = PMMG2D_locatePoint( initialMesh, &list_triangles2[idx][idy][0], ppt, barycoord,
                                    &ifoundTria1, &ifoundEdge1, &ifoundVertex1 );

          if( mesh->info.imprim > PMMG2D_VERB_ITWAVES )
            PMMG2D_locatePoint_errorCheck( mesh,ip,ier,myrank );

          // Interpolate point metrics
          if( ismet ) {
            if( ifoundVertex1 != PMMG2D_UNSET ) {
              ier = PMMG2D_copyMetrics( mesh,met,initialMesh,initialMet,ip,
                                        initialMesh->tria[ifoundTria1].v[ifoundVertex1] );
            } else if( ifoundEdge1 != PMMG2D_UNSET ) {
              ier = PMMG2D_interp2bar( mesh,met,initialMet,&initialMesh->tria[ifoundTria1],
                                       ip,ifoundEdge1,barycoord );
            } else if( !nlayers ) {
              ier = PMMG2D_interp3bar(mesh,met,initialMet,&initialMesh->tria[ifoundTria1],ip,
                                      barycoord);
            } else {
              ier = PMMG2D_extended_interpNbar(initialMesh,met,initialMet,nlayers,ppt->c[0],ppt->c[1],ifoundTria1,ip);
            }
          }
        }
        else {

          // Find the cell containing the point
          idx = (int)((ppt->c[0] - xmin) / ((xmax - xmin)/GRID_SIZE));
          if (idx == -1) idx = 0; // Possible if the point is on the domaine boundary
          if (idx == GRID_SIZE) idx = GRID_SIZE-1;
          idy = (int)((ppt->c[1] - ymin) / ((ymax - ymin)/GRID_SIZE));
          if (idy == -1) idy = 0;
          if (idy == GRID_SIZE) idy = GRID_SIZE-1;

          // Locate point in the old mesh 
          ier = PMMG2D_locatePoint( oldMesh, &list_triangles[idx][idy][0], ppt, barycoord,
                                    &ifoundTria2, &ifoundEdge2, &ifoundVertex2 );

          if( mesh->info.imprim > PMMG2D_VERB_ITWAVES )
            PMMG2D_locatePoint_errorCheck( mesh,ip,ier,myrank );
  
          // Interpolate point metrics
          if( ismet ) {
            if( ifoundVertex2 != PMMG2D_UNSET ) {
              ier = PMMG2D_copyMetrics( mesh,met,oldMesh,oldMet,ip,
                                        oldMesh->tria[ifoundTria2].v[ifoundVertex2] );
            } else if( ifoundEdge2 != PMMG2D_UNSET ) {
              ier = PMMG2D_interp2bar( mesh,met,oldMet,&oldMesh->tria[ifoundTria2],
                                       ip,ifoundEdge2,barycoord );
            } else if ( !nlayers ) {
              ier = PMMG2D_interp3bar(mesh,met,oldMet,&oldMesh->tria[ifoundTria2],ip,
                                      barycoord);
            } else {
              ier = PMMG2D_extended_interpNbar(oldMesh,met,oldMet,nlayers,ppt->c[0],ppt->c[1],ifoundTria2,ip);
            }
          }
        }

        // Flag point as interpolated
        ppt->flag = mesh->base;

      }
    }
  }

  // Free the pointers
  for (i = 0; i < GRID_SIZE; i++) {
    for (k = 0; k < GRID_SIZE; k++) {
      free(list_triangles[i][k]);
    }
    free(list_triangles[i]);
  }
  free(list_triangles);

  if (!optim_interp) return 1;

  for (i = 0; i < GRID_SIZE; i++) {
    for (k = 0; k < GRID_SIZE; k++) {
      free(list_triangles2[i][k]);
    }
    free(list_triangles2[i]);
  }
  free(list_triangles2);

  return 1;

}

/**
 * \param parmesh pointer to the parmesh structure.
 * \param field array to interpolate
 * \return 0 if fail, 1 if success
 *
 * Interpolate a field from background to current meshes
 * and store the inteprolation in the normal
 *
 */
int PMMG2D_interpFields( PMMG2D_pParMesh parmesh, double* field ) {
  MMG5_pTria pt;
  MMG5_pPoint ppt;
  MMG5_pMesh mesh = parmesh->mesh;
  MMG5_pMesh oldMesh = parmesh->old_mesh;

  PMMG2D_barycoord barycoord[3];
  int         ifoundTria;
  int         ifoundEdge;
  int         ifoundVertex;
  int         ip, ie, iloc, idx, idy, i, k;
  int         ier;
  int***      list_triangles;
  int         GRID_SIZE = sqrt(mesh->nt);
  double      xmin, xmax, ymin, ymax;

  find_bounding_box(oldMesh, &xmin, &xmax, &ymin, &ymax);
  list_triangles = grid_size_triangles(oldMesh, xmin, ymin, xmax, ymax, GRID_SIZE);

  ifoundTria = 1;

  // Loop on new triangles, and localize their vertices in the old mesh
  for( ie = 1; ie <= mesh->nt; ie++ ) {
    pt = &mesh->tria[ie];

    if( !MG_EOK(pt) ) continue;

    for( iloc = 0; iloc < 3; iloc++ ) {
      ip = pt->v[iloc];
      ppt = &mesh->point[ip];

      if( !MG_VOK(ppt) ) continue;

      // Find the cell containing the point
      idx = (int)((ppt->c[0] - xmin) / ((xmax - xmin)/GRID_SIZE));
      if (idx == -1) idx = 0; // Possible if the point is on the domaine boundary
      if (idx == GRID_SIZE) idx = GRID_SIZE-1;
      idy = (int)((ppt->c[1] - ymin) / ((ymax - ymin)/GRID_SIZE));
      if (idy == -1) idy = 0;
      if (idy == GRID_SIZE) idy = GRID_SIZE-1;

      // Locate point in the old mesh 
      ier = PMMG2D_locatePoint( oldMesh, &list_triangles[idx][idy][0], ppt, barycoord,
                                &ifoundTria, &ifoundEdge, &ifoundVertex );

      // Interpolate point metrics
      if( ifoundVertex != PMMG2D_UNSET ) {
        ppt->n[0] = field[oldMesh->tria[ifoundTria].v[ifoundVertex]-1];
        ppt->n[1] = field[oldMesh->np + oldMesh->tria[ifoundTria].v[ifoundVertex]-1];
      }
      else if( ifoundEdge != PMMG2D_UNSET ) {

        int    i0,i1;
        double phi[3];
        PMMG2D_barycoord_get( phi, barycoord, 3 );

        i0 = MMG5_inxt2[ifoundEdge];
        i1 = MMG5_iprv2[ifoundEdge];

        // Linear interpolation of the squared size
        ppt->n[0] = phi[i0]*field[oldMesh->tria[ifoundTria].v[i0]-1] +
                   phi[i1]*field[oldMesh->tria[ifoundTria].v[i1]-1];
        ppt->n[1] = phi[i0]*field[oldMesh->np + oldMesh->tria[ifoundTria].v[i0]-1] +
                   phi[i1]*field[oldMesh->np + oldMesh->tria[ifoundTria].v[i1]-1];
      }
      else {
        double phi[3];
        int i,j;

        PMMG2D_barycoord_get( phi, barycoord, 3 );

        // Linear interpolation of the squared size
        ppt->n[0] = 0.;
        ppt->n[1] = 0.;
        for( i = 0; i < 3; i++ ) {
          // Barycentric coordinates could be permuted
          ppt->n[0] += phi[i]*field[oldMesh->tria[ifoundTria].v[i]-1];
          ppt->n[1] += phi[i]*field[oldMesh->tria[ifoundTria].v[i]-1 + oldMesh->np];
        }
      }
    }
  }

  // Free the pointers
  for (i = 0; i < GRID_SIZE; i++) {
    for (k = 0; k < GRID_SIZE; k++) {
      free(list_triangles[i][k]);
    }
    free(list_triangles[i]);
  }
  free(list_triangles);

  return 1;
}


/**
 * \param parmesh pointer to the parmesh structure.
 *
 * \return 0 if fail, 1 if success
 *
 *  Interpolate metrics from background to current meshes.
 *  Do nothing if no metrics is provided (info.inputMet == 0), otherwise:
 *  - if the metrics is constant, recompute it;
 *  - else, interpolate the non-constant metrics.
 *
 */
int PMMG2D_interpMetrics( PMMG2D_pParMesh parmesh, double** polygon, int size_polygon ) {
  int ier = 1;

  if( !PMMG2D_interpMetrics_mesh( parmesh->mesh, parmesh->old_mesh,
                                  parmesh->initial_mesh,
                                  parmesh->met, parmesh->old_met,
                                  parmesh->initial_met,
                                  polygon, size_polygon,
                                  parmesh->info.inputMet,
                                  parmesh->info.optim_interp,
                                  parmesh->info.interp_layers,
                                  parmesh->myrank) ) ier = 0;

  return ier;
}

/**
 * \param mesh pointer to the current mesh.
 * \param size pointer to the size of the output contour.
 *
 * \return the contour of the partition
 *
 * Compute the polygon which is the contour of the partition.
 *
 */
double** build_partition_contour(MMG5_pMesh mesh, int *size) {
  int k, i, cur, next;
  int first = 0;
  int* list_points = (int*) malloc((mesh->np+1) * sizeof(int));
  MMG5_pTria pt;

  for (k = 0; k <= mesh->np; k++) list_points[k] = 0;

  // Find the boundary and parallel interface nodes
  for (k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if( !MG_EOK(pt) ) continue;

    for (i = 0; i < 3; i++) {
      if ((pt->tag[i] & MG_BDY) || (pt->tag[i] & MG_PARBDY)) list_points[pt->v[(i+1)%3]] = pt->v[(i+2)%3];
      if (!first && (pt->tag[i] & MG_PARBDY)) first = pt->v[(i+1)%3];
    }
  }

  k = 0;
  cur = first;
  next = 0;
  while (first != next) {
    next = list_points[cur];
    cur = next;
    k++;
  }
  *size = k;

  double **contour = (double**)malloc(*size * sizeof(double *));
  
  cur = first;
  for (k = 0; k < *size; k++) {
    contour[k] = (double*)malloc(2*sizeof(double));
    contour[k][0] = mesh->point[cur].c[0];
    contour[k][1] = mesh->point[cur].c[1];
    cur = list_points[cur];
  }

  free(list_points);

  return contour;
}

static int min(int* list) {
  int res = list[0];
  for (int i = 1; i <= 5; i++) if (list[i] < res) res = list[i];
  return res;
}

static int max(int* list) {
  int res = list[0];
  for (int i = 1; i <= 5; i++) if (list[i] > res) res = list[i];
  return res;
}

/**
 * \param mesh pointer to the current mesh.
 * \param minX x-coordinate of the leftmost point of the mesh.
 * \param minY y-coordinate of the lowest point of the mesh.
 * \param maxX x-coordinate of the rightmost point of the mesh.
 * \param maxY y-coordinate of the highest point of the mesh.
 * \param GRID_SIZE the number of cell in each direction.
 *
 * \return the list of triangles in each cell
 *
 * Divide the domain as a Cartesian grid and store the triangles in each cell of this grid.
 *
 */
int*** grid_size_triangles(MMG5_pMesh mesh, double minX, double minY, double maxX, double maxY, int GRID_SIZE) {
  MMG5_pTria pt;
  int i, k;
  int idx, idy, maxx, maxy, minx, miny;
  int list_idx[6], list_idy[6];

  double cellWidth = (maxX - minX) / GRID_SIZE;
  double cellHeight = (maxY - minY) / GRID_SIZE;

  int ***list_triangles = (int***)malloc(GRID_SIZE * sizeof(int**));
  for (i = 0; i < GRID_SIZE; i++) {
    list_triangles[i] = (int**)malloc(GRID_SIZE * sizeof(int*));
    for (k = 0; k < GRID_SIZE; k++) {
      list_triangles[i][k] = (int*)malloc(sizeof(int));
      list_triangles[i][k][0] = 0;
    }
  }

  for (k = 1; k <= mesh->nt; k++) {

    pt = &mesh->tria[k];
    
    if ( !MG_EOK(pt) )  continue;

    // List of cells containing at least a part of the triangle
    for (i = 0; i < 3; i++) {
      list_idx[i] = (int)((mesh->point[pt->v[i]].c[0] - minX) / cellWidth + PMMG2D_EPS);
      list_idy[i] = (int)((mesh->point[pt->v[i]].c[1] - minY) / cellHeight + PMMG2D_EPS);
      list_idx[3+i] = (int)((mesh->point[pt->v[i]].c[0] - minX) / cellWidth - PMMG2D_EPS);
      list_idy[3+i] = (int)((mesh->point[pt->v[i]].c[1] - minY) / cellHeight - PMMG2D_EPS);
    }

    minx = min(list_idx);
    if (minx < 0) minx = 0;
    maxx = max(list_idx);
    if (maxx >= GRID_SIZE) maxx = GRID_SIZE-1;
    miny = min(list_idy);
    if (miny < 0) miny = 0;
    maxy = max(list_idy);
    if (maxy >= GRID_SIZE) maxy = GRID_SIZE-1;

    for (idx = minx; idx <= maxx; idx++) {
      for (idy = miny; idy <= maxy; idy++) {
        list_triangles[idx][idy][0]++;
        list_triangles[idx][idy] = (int*)realloc(list_triangles[idx][idy], (list_triangles[idx][idy][0]+1)*sizeof(int));
        list_triangles[idx][idy][list_triangles[idx][idy][0]] = k;
      }
    }
  }

  return list_triangles;
}
