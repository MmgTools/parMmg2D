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
 * \file distributemesh_pmmg2d.c
 * \brief Distribute the mesh over the processors.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */
#include "parmmg2d.h"
#include "mpitypes_pmmg2d.h"
#include "metis_pmmg2d.h"

static int compare(const void *a, const void *b) {
  return (*(int*)a - *(int*)b);
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 * \param point_structure contains the points of the mesh.
 * \param m contains the metric values over the mesh.
 *
 * Renumber the vertices in the local meshes.
 */
void PMMG2D_renumber_vertices( PMMG2D_pParMesh parmesh, 
                               MMG5_pPoint point_structure,
                               double* m)
{
  MMG5_pTria pt;
  int* list_points = (int*) malloc(3*parmesh->mesh->nt * sizeof(int));
  int* conversion_table_vertices = (int*) malloc((parmesh->mesh->np+1) * sizeof(int));
  int* inverse_conversion_table_vertices = (int*) malloc((parmesh->mesh->np+1) * sizeof(int));
  int j = 0;
  int l = 2;
  int ier;

  // 1) Create the list of points taking part in the triangles
  for (int k = 1; k <= parmesh->mesh->nt; k++) {
    pt = &parmesh->mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (int i = 0; i < 3; i++) {
      list_points[j] = pt->v[i];
      j++;
    }
  }

  // 2) Sort this list in ascending order
  qsort(list_points, j, sizeof(int), compare);

  // 3) Associate local indices to global indices
  conversion_table_vertices[list_points[0]] = 1;
  inverse_conversion_table_vertices[1] = list_points[0];
  for (int k = 1; k < j; k++) {
    if (list_points[k] == list_points[k-1]) continue;
    conversion_table_vertices[list_points[k]] = l;
    inverse_conversion_table_vertices[l] = list_points[k];
    l++;
  }

  // 4) Renumber vertices in triangles and points
  for (int k = 1; k <= parmesh->mesh->nt; k++) {
    pt = &parmesh->mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    for (int i = 0; i < 3; i++) pt->v[i] = conversion_table_vertices[pt->v[i]];
  }

  parmesh->mesh->np = l-1;
  parmesh->mesh->npi = parmesh->mesh->np;
  parmesh->mesh->npnil = 0;
  parmesh->met->np = parmesh->mesh->np;
  parmesh->met->npi = parmesh->mesh->np;
  parmesh->mesh->npmax = parmesh->mesh->np;
  parmesh->met->npmax = parmesh->met->np;

  if (!parmesh->myrank) {
    PMMG2D_DEL_MEM(parmesh->mesh,parmesh->mesh->point, MMG5_Point, "deallocate points");
  }
  PMMG2D_CALLOC(parmesh->mesh,parmesh->mesh->point,parmesh->mesh->np+1,MMG5_Point,"parallel points", ier = 0);
  PMMG2D_CALLOC(parmesh->mesh,parmesh->met->m,parmesh->met->size*(parmesh->met->np+1),
                double, "metric parallel", ier = 6);

  for (int k = 1; k <= parmesh->mesh->np; k++) {
    parmesh->mesh->point[k] = point_structure[inverse_conversion_table_vertices[k]];
    if (parmesh->met->size == 1) { // isotropic
      parmesh->met->m[k] = m[inverse_conversion_table_vertices[k]];
    }
    else {
      for (int i = 0; i < 3; i++) { // anisotropic
        parmesh->met->m[3*k+i] = m[3*inverse_conversion_table_vertices[k]+i];
      }
    }
  }

  free(list_points);
  free(conversion_table_vertices);
  free(inverse_conversion_table_vertices);
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Broadcast the options to the other processes
 *
 */
int PMMG2D_Bcast_info( PMMG2D_pParMesh parmesh )
{
  int ier;

  MMG5_pMesh mesh = parmesh->mesh;

  MPI_CHECK( MPI_Bcast( &mesh->memMax,         1, MPI_LONG_LONG, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.dhd,       1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hmin,      1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hmax,      1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hsiz,      1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hgrad,     1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hgradreq,  1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.hausd,     1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);

  MPI_CHECK( MPI_Bcast( &mesh->info.delta,     1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.min,       3, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);

  MPI_CHECK( MPI_Bcast( &mesh->info.ls,        1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);

  MPI_CHECK( MPI_Bcast( &mesh->info.npar,      1, MPI_INT, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.opnbdy,    1, MPI_INT, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.renum,     1, MPI_INT, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.PROctree,  1, MPI_INT, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.nmat,      1, MPI_INT, parmesh->info.root, parmesh->comm ), ier=6);

  MPI_CHECK( MPI_Bcast( &mesh->info.nreg,      1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.nosizreq,  1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.imprim,    1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.ddebug,    1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.iso,       1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.lag,       1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.parTyp,    1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.optim,     1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.noinsert,  1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.noswap,    1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.nomove,    1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
  MPI_CHECK( MPI_Bcast( &mesh->info.nosurf,    1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);

  if ( mesh->info.npar ) {

    if( parmesh->myrank != parmesh->info.root )
      MMG5_SAFE_CALLOC(mesh->info.par,mesh->info.npar,MMG5_Par, ier = 0);

    if ( ier ) { 
      for ( int k = 0; k < mesh->info.npar; ++k ) {
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].hmin,    1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].hmax,    1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].hausd,   1, MPI_DOUBLE, parmesh->info.root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].ref,     1, MPI_INT, parmesh->info.root, parmesh->comm ), ier=6);
        MPI_CHECK( MPI_Bcast( &mesh->info.par[k].elt,     1, MPI_CHAR, parmesh->info.root, parmesh->comm ), ier=6);
      }
    }
  }

  mesh->base = 0;

  return ier;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 *
 * \return 0 (on all procs) if fail, 1 otherwise
 *
 * Send the local meshes to the corresponding procs.
 */
int PMMG2D_exchange_mesh( PMMG2D_pParMesh parmesh )
{
  // Note: the edges are destroyed by MMG2D, so they cannot be exchanged here (so edge=false)
  int ier,k,n,i,na_tot;
  int list_index[parmesh->nprocs];
  int edge = 0;
  int met = 0;
  if (parmesh->met && parmesh->met->m) met = 1;
  if (parmesh->mesh->edge && parmesh->mesh->na) edge = 1;

  MPI_CHECK( MPI_Bcast( &met, 1, MPI_INT, parmesh->info.root, parmesh->comm ), ier = 0);
  MPI_CHECK( MPI_Bcast( &edge, 1, MPI_INT, parmesh->info.root, parmesh->comm ), ier = 0);

  MMG5_pTria pt_global;
  MMG5_pTria* pt = NULL;
  MMG5_pPoint* pp = NULL;
  MMG5_pEdge* ped = NULL;

  MPI_Datatype   MPI_POINT, MPI_TRIANGLE, MPI_EDGE;
  MPI_Status     status;

  // First build the exchange array on rank 0
  if (parmesh->myrank == parmesh->info.root ) {

    pt = (MMG5_pTria*) malloc(parmesh->nprocs * sizeof(MMG5_pTria));
    for (k = 0; k < parmesh->nprocs; k++) pt[k] = (MMG5_pTria) malloc(parmesh->mesh->nt * sizeof(MMG5_Tria));

    pp = (MMG5_pPoint*) malloc(parmesh->nprocs * sizeof(MMG5_pPoint));
    for (k = 0; k < parmesh->nprocs; k++) pp[k] = (MMG5_pPoint) malloc((parmesh->mesh->np+1) * sizeof(MMG5_Point));

    if (edge) {
      ped = (MMG5_pEdge*) malloc(parmesh->nprocs * sizeof(MMG5_pEdge));
      for (k = 0; k < parmesh->nprocs; k++) ped[k] = (MMG5_pEdge) malloc(parmesh->mesh->na * sizeof(MMG5_Edge));
    }

    for (k = 0; k < parmesh->nprocs; k++) list_index[k] = 0;

    // Dispatch triangles, vertices and edges in arrays relative
    // to each process
    for (k = 1; k <= parmesh->mesh->nt; k++) {
      pt_global = &parmesh->mesh->tria[k];
      if ( !MG_EOK(pt_global) ) continue;

      n = pt_global->base;

      pt[n][list_index[n]] = parmesh->mesh->tria[k];

      // Add triangle vertices and edges in point and edge temporary array respectively
      for ( i = 0; i < 3 ; ++i ) {
        pp[n][pt_global->v[i]] = parmesh->mesh->point[ pt_global->v[i] ];
        if (edge) ped[n][pt_global->edg[i]] = parmesh->mesh->edge[ pt_global->edg[i] ];
      }

      // Increment the triangle index
      list_index[n]++;
    }
  }

  // Create specific types to MPI exchanges
  PMMG2D_create_MPI_Tria( &MPI_TRIANGLE );  
  if (edge) PMMG2D_create_MPI_Edge( &MPI_EDGE );
  PMMG2D_create_MPI_Point( &MPI_POINT );

  // Send triangles, edges and vertices to other processes
  // Initially, consider the global number of points locally
  MPI_CHECK( MPI_Bcast( &parmesh->mesh->np, 1, MPI_INT, parmesh->info.root, parmesh->comm ), ier=6);
  if (met) {
    MPI_CHECK( MPI_Bcast( &parmesh->met->size, 1, MPI_INT, parmesh->info.root, parmesh->comm ), ier=6);
    if ( parmesh->met->size == 1 ) {
      parmesh->met->type = MMG5_Scalar;
    }
    else if ( parmesh->met->size == 3 ) {
      parmesh->met->type = MMG5_Tensor;
    }
  }

  if (parmesh->myrank == parmesh->info.root) {

    for ( n = 0; n < parmesh->nprocs; n++ ) {

      if (n == parmesh->info.root) continue;

      // Triangles
      MPI_CHECK( MPI_Send(&list_index[n], 1, MPI_INT, n, MPI_ANALYS_TAG+n,
                           parmesh->comm), ier = 0 );

      MPI_CHECK( MPI_Send(pt[n], list_index[n], MPI_TRIANGLE, n, MPI_ANALYS_TAG + parmesh->nprocs + 1 + n,
                           parmesh->comm), ier = 0 );

      // Edges
      if (edge) {
        MPI_CHECK( MPI_Send(&parmesh->mesh->na, 1, MPI_INT, n, MPI_ANALYS_TAG + 2*parmesh->nprocs + 1 + n,
                             parmesh->comm), ier = 0 );
        MPI_CHECK( MPI_Send(ped[n], parmesh->mesh->na, MPI_EDGE, n, MPI_ANALYS_TAG + 3*parmesh->nprocs + 1 + n,
                             parmesh->comm), ier = 0 );
      }

      // Points
      MPI_CHECK( MPI_Send(pp[n], parmesh->mesh->np+1, MPI_POINT, n, MPI_ANALYS_TAG + 4*parmesh->nprocs + 1 + n,
                           parmesh->comm), ier = 0 );

      // Metric
      if (met) MPI_CHECK( MPI_Send(&parmesh->met->m[0], parmesh->met->size*(parmesh->mesh->np+1), MPI_DOUBLE, n, 
                                   MPI_ANALYS_TAG + 5*parmesh->nprocs + 1 + n, parmesh->comm), ier = 0 );
    }

    // Deal with rank 0 information
    PMMG2D_DEL_MEM(parmesh->mesh,parmesh->mesh->tria, MMG5_Tria, "deallocate triangles");

    parmesh->mesh->nt = list_index[parmesh->info.root];
    parmesh->mesh->nenil = 0;
    parmesh->mesh->ntmax = parmesh->mesh->nti = parmesh->mesh->nt;
    parmesh->mesh->xtmax = parmesh->mesh->ntmax;

    PMMG2D_CALLOC(parmesh->mesh,parmesh->mesh->tria,parmesh->mesh->nt+1,MMG5_Tria,"parallel triangles", return 0);

    MMG5_pPoint pp_local = (MMG5_pPoint) malloc((parmesh->mesh->np+1) * sizeof(MMG5_Point));
    for ( n = 0; n < parmesh->mesh->nt; n++) parmesh->mesh->tria[n+1] = pt[parmesh->info.root][n];

    for ( n = 1; n < parmesh->mesh->np+1; n++) pp_local[n] = pp[parmesh->info.root][n];

    for (k = 0; k < parmesh->nprocs; k++) {
      free(pt[k]);
      free(pp[k]);
      if (edge) free(ped[k]);
    }
    free(pt);
    free(pp);
    if (edge) free(ped);

    // Renumber the vertices with local indices
    PMMG2D_renumber_vertices(parmesh, pp_local, parmesh->met->m);
    free(pp_local);
  }
  else {
    // Triangles
    MPI_CHECK( MPI_Recv(&parmesh->mesh->nt, 1, MPI_INT, parmesh->info.root, 
                         MPI_ANALYS_TAG + parmesh->myrank,
                         parmesh->comm, &status), ier = 0 );

    parmesh->mesh->ntmax = parmesh->mesh->nti = parmesh->mesh->nt;
    parmesh->mesh->nenil = 0;
    parmesh->mesh->xtmax = parmesh->mesh->ntmax;

    PMMG2D_CALLOC(parmesh->mesh,parmesh->mesh->tria,parmesh->mesh->nt+1,MMG5_Tria,"parallel triangles", return 0);

    MMG5_pTria pt_local = (MMG5_pTria) malloc((parmesh->mesh->nt+1) * sizeof(MMG5_Tria));
    MPI_CHECK( MPI_Recv(pt_local, parmesh->mesh->nt+1, MPI_TRIANGLE, parmesh->info.root, 
                        MPI_ANALYS_TAG + parmesh->myrank + parmesh->nprocs + 1,
                        parmesh->comm, &status), ier = 0 );

    for (int k = 1; k <= parmesh->mesh->nt; k++) parmesh->mesh->tria[k] = pt_local[k-1];

    free(pt_local);

    // Edges
    if (edge) {
      MPI_CHECK( MPI_Recv(&na_tot, 1, MPI_INT, parmesh->info.root,
                           MPI_ANALYS_TAG + 2*parmesh->nprocs + 1 + parmesh->myrank,
                           parmesh->comm, &status), ier = 0 );
      MMG5_pEdge ped_local = (MMG5_pEdge) malloc(na_tot * sizeof(MMG5_Edge));
      MPI_CHECK( MPI_Recv(ped_local, na_tot, MPI_EDGE, parmesh->info.root,
                          MPI_ANALYS_TAG + 3*parmesh->nprocs + 1 + parmesh->myrank,
                          parmesh->comm, &status), ier = 0 );
      free(ped_local);
    }

    // Points 
    MMG5_pPoint pp_local = (MMG5_pPoint) malloc((parmesh->mesh->np+1) * sizeof(MMG5_Point));
    MPI_CHECK( MPI_Recv(pp_local, parmesh->mesh->np+1, MPI_POINT, parmesh->info.root,
                        MPI_ANALYS_TAG + 4*parmesh->nprocs + 1 + parmesh->myrank,
                        parmesh->comm, &status), ier = 0 );

    // Metric
    double* met_global = (double*) malloc(parmesh->met->size*(parmesh->mesh->np+1) * sizeof(double));
    if (met) MPI_CHECK( MPI_Recv(met_global, parmesh->met->size*(parmesh->mesh->np+1), MPI_DOUBLE, parmesh->info.root,
                                 MPI_ANALYS_TAG + 5*parmesh->nprocs + 1 + parmesh->myrank,
                                 parmesh->comm, &status), ier = 0 );

    // Renumber the vertices with local indices
    PMMG2D_renumber_vertices(parmesh, pp_local, met_global);
    free(pp_local);
    free(met_global);
  }

  ier = PMMG2D_Bcast_info(parmesh);

  PMMG2D_Free_MPI_meshDatatype(&MPI_POINT, &MPI_EDGE, &MPI_TRIANGLE, edge);

  return 1;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 *
 * \return 0 (on all procs) if fail, 1 otherwise
 *
 * Send the local meshes from root to the corresponding procs.
 */
int PMMG2D_exchange_from_root( PMMG2D_pParMesh parmesh, int size, MMG5_pTria* pt_global, MMG5_pPoint* pp_global, double** mm_global, int* npt_global, int* npp_global )
{
  int n, npt = 0, npp = 0;
  int ier = 1;
  int size_met = parmesh->met->size;

  MPI_Datatype   MPI_POINT, MPI_TRIANGLE;
  PMMG2D_create_MPI_Tria( &MPI_TRIANGLE );
  PMMG2D_create_MPI_Point( &MPI_POINT );
  MPI_Request send_request, recv_request;
  MPI_Status     status;
  
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;
  MMG5_pTria pt;

  *npt_global = 0;
  *npp_global = 0;
  
  for ( n = 1; n < parmesh->nprocs; n++ ) {

    if (parmesh->myrank != n && parmesh->myrank != parmesh->info.root) continue;

    int length = 0;
    for (int k = 1; k <= mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;
      if (pt->base != n) continue;
      length++;
    }

    int index_pt = 0;
    int index_pp = 0;

    MMG5_pTria list_pt = NULL;
    MMG5_pPoint list_pp = NULL;
    double* list_mm = NULL;
    if (parmesh->myrank == parmesh->info.root) {
      list_pt = (MMG5_pTria) malloc(length * sizeof(MMG5_Tria));
      list_pp = (MMG5_pPoint) malloc(3*length * sizeof(MMG5_Point));
      list_mm = (double*) malloc(parmesh->met->size*3*length * sizeof(double));
    }

    // 1) Find the triangles to exchange
    for (int k = 1; k <= mesh->nt; k++) {

      if (parmesh->myrank != parmesh->info.root) continue;

      pt = &mesh->tria[k];

      if ( !MG_EOK(pt) ) continue;

      if (pt->base != n) continue;

      list_pt[index_pt++] = mesh->tria[k];

      for (int i = 0; i <= 2; i++) {
        list_pp[index_pp] = mesh->point[pt->v[i]];
        if (parmesh->met->size == 1) { // isotropic
          list_mm[index_pp] = parmesh->met->m[pt->v[i]];
        }
        else { // anisotropic
          list_mm[3*index_pp] = parmesh->met->m[3*pt->v[i]];
          list_mm[3*index_pp+1] = parmesh->met->m[3*pt->v[i]+1];
          list_mm[3*index_pp+2] = parmesh->met->m[3*pt->v[i]+2];
        }
        index_pp++;
      }
    }

    // Triangles
    if (parmesh->myrank == parmesh->info.root) MPI_CHECK( MPI_Send(&index_pt, 1, MPI_INT, n, MPI_ANALYS_TAG + n, parmesh->comm), ier = 0 );

    if (parmesh->myrank == n) MPI_CHECK( MPI_Recv(&npt, 1, MPI_INT, parmesh->info.root, MPI_ANALYS_TAG + parmesh->myrank,
                                                  parmesh->comm, &status), ier = 0 );

    if (parmesh->myrank == parmesh->info.root) MPI_CHECK( MPI_Isend(list_pt, index_pt, MPI_TRIANGLE, n,
                                                          MPI_ANALYS_TAG + n + parmesh->nprocs, parmesh->comm, &send_request), ier = 0 );
    MMG5_pTria pt_local = NULL; 
    if (parmesh->myrank == n) {
      *npt_global += npt;
      pt_local = (MMG5_pTria) malloc(npt * sizeof(MMG5_Tria));
      MPI_CHECK( MPI_Irecv(pt_local, npt, MPI_TRIANGLE, parmesh->info.root,
                           MPI_ANALYS_TAG + parmesh->myrank + parmesh->nprocs,
                           parmesh->comm, &recv_request), ier = 0 );
    }

    if (parmesh->myrank == parmesh->info.root) MPI_Wait(&send_request, &status);
    if (parmesh->myrank == n) MPI_Wait(&recv_request, &status);

    if (parmesh->myrank == n) {
      PMMG2D_REALLOC(mesh, *pt_global, *npt_global, *npt_global-npt,
                     MMG5_Tria,"Front advancing exchange triangle",return 0);
      memcpy(*pt_global + *npt_global - npt, pt_local, npt*sizeof(MMG5_Tria));
      free(pt_local);
    }

    // Points
    if (parmesh->myrank == parmesh->info.root) MPI_CHECK( MPI_Send(&index_pp, 1, MPI_INT, n, MPI_ANALYS_TAG + n + 2*parmesh->nprocs,
                                                                   parmesh->comm), ier = 0 );

    if (parmesh->myrank == n) MPI_CHECK( MPI_Recv(&npp, 1, MPI_INT, parmesh->info.root, MPI_ANALYS_TAG + parmesh->myrank + 2*parmesh->nprocs,
                                                  parmesh->comm, &status), ier = 0 );

    if (parmesh->myrank == parmesh->info.root) MPI_CHECK( MPI_Isend(list_pp, index_pp, MPI_POINT, n,
                                                          MPI_ANALYS_TAG + n + 3*parmesh->nprocs, parmesh->comm, &send_request), ier = 0 );
    MMG5_pPoint pp_local = NULL;
    if (parmesh->myrank == n) {
      *npp_global += npp;
      pp_local = (MMG5_pPoint) malloc(npp*sizeof(MMG5_Point));
      MPI_CHECK( MPI_Irecv(pp_local, npp, MPI_POINT, parmesh->info.root,
                           MPI_ANALYS_TAG + parmesh->myrank + 3*parmesh->nprocs,
                           parmesh->comm, &recv_request), ier = 0 );
    }

    if (parmesh->myrank == parmesh->info.root) MPI_Wait(&send_request, &status);
    if (parmesh->myrank == n) MPI_Wait(&recv_request, &status);

    if (parmesh->myrank == n) {
      PMMG2D_REALLOC(mesh, *pp_global, *npp_global, *npp_global - npp,
                     MMG5_Point,"Front advancing exchange point",return 0);

      memcpy(*pp_global + *npp_global - npp, pp_local, npp*sizeof(MMG5_Point));
      free(pp_local);
    }

    // Metrics
    if (parmesh->myrank == parmesh->info.root) MPI_CHECK( MPI_Isend(list_mm, size_met*index_pp, MPI_DOUBLE, n,
                                                                    MPI_ANALYS_TAG + n + 4*parmesh->nprocs, parmesh->comm, &send_request), ier = 0 );

    MMG5_pSol mm_local = NULL;
    if (parmesh->myrank == n) {
      mm_local = (MMG5_pSol) malloc(size_met*npp * sizeof(MMG5_Sol));
      MPI_CHECK( MPI_Irecv(mm_local, size_met*npp, MPI_DOUBLE, parmesh->info.root,
                           MPI_ANALYS_TAG + parmesh->myrank + 4*parmesh->nprocs,
                           parmesh->comm, &recv_request), ier = 0 );
    }

    if (parmesh->myrank == parmesh->info.root) MPI_Wait(&send_request, &status);
    if (parmesh->myrank == n) MPI_Wait(&recv_request, &status);

    if (parmesh->myrank == n) {
      PMMG2D_REALLOC(mesh, *mm_global, size_met* *npp_global, size_met* (*npp_global - npp),
                     double,"Front advancing exchange metric",return 0);

      memcpy(*mm_global + size_met*(*npp_global - npp), mm_local, size_met*npp*sizeof(double));
      free(mm_local);
    }
  
    if (parmesh->myrank == parmesh->info.root) {
      free(list_pt);
      free(list_pp);
      free(list_mm);
    }

  }
  
  PMMG2D_Free_MPI_meshDatatype(&MPI_POINT, &MPI_TRIANGLE, &MPI_TRIANGLE, false);

  return ier;
}

// Write vtk to check the mpi exchange
static void check_vtk_file(PMMG2D_pParMesh parmesh)
{
  MMG5_pTria pt;
  char name_file[20];
  sprintf(name_file, "mesh_%d.vtk", parmesh->myrank);
 
  // Create the list of parallel boundary points
  PMMG2D_fill_interface_nodes_list( parmesh );

  FILE *fp = fopen(name_file, "w");
  fprintf(fp, "# vtk DataFile Version 3.0\nMesh with metrics at points\nASCII\nDATASET UNSTRUCTURED_GRID\n");
  fprintf(fp, "POINTS %d double\n", parmesh->mesh->np);
  
  // Vertices
  for (int n = 1; n < parmesh->mesh->np+1; n++) { 
    fprintf(fp, "%f %f 0.\n", parmesh->mesh->point[n].c[0], parmesh->mesh->point[n].c[1]);
  }

  fprintf(fp, "\nCELLS %d %d\n", parmesh->mesh->nt, 4*parmesh->mesh->nt);

  // Triangles
  for (int k = 1; k < parmesh->mesh->nt+1; k++) {
    pt = &parmesh->mesh->tria[k];
    fprintf(fp, "3 %d %d %d \n",pt->v[0]-1, pt->v[1]-1, pt->v[2]-1);
  }

  fprintf(fp, "\nCELL_TYPES %d\n", parmesh->mesh->nt);
  for (int k = 1; k < parmesh->mesh->nt+1; k++) fprintf(fp, "5\n"); 

  if (!parmesh->met && !parmesh->met->m) {
    fclose(fp);
    return;
  }

  fprintf(fp, "\nPOINT_DATA %d\n", parmesh->mesh->np);

  // Parallel boundary points
  int is_boundary_nodes[parmesh->mesh->np];
  for (int n = 0; n < parmesh->mesh->np; n++) is_boundary_nodes[n] = 0;
  for (int n = 0; n < parmesh->nip; n++) 
    is_boundary_nodes[(&parmesh->list_interface_nodes[n])->index[0]-1] = (&parmesh->list_interface_nodes[n])->proc[1]+1;

  // SCALARS fail when reading with paraview 5.10 for unknown reason
  fprintf(fp, "VECTORS Boundary int \n");
  for (int n = 0; n < parmesh->mesh->np; n++)
    fprintf(fp, "%d %d %d\n", is_boundary_nodes[n], 0, 0);

  fprintf(fp, "\n");

  // Metrics
  fprintf(fp, "VECTORS Metric double \n");
  if (parmesh->met->size == 1) {
    for (int n = 1; n < parmesh->mesh->np+1; n++) {
      fprintf(fp, "%f %f %f\n", parmesh->met->m[n], 0., 0.);
    }
  }
  else {
    for (int n = 1; n < parmesh->mesh->np+1; n++) { 
      fprintf(fp, "%f %f %f\n", parmesh->met->m[3*n], parmesh->met->m[3*n+1], parmesh->met->m[3*n+2]);
    }
  }

  fclose(fp);

}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 *
 * \return 0 (on all procs) if fail, 1 otherwise
 *
 * Send the local meshes to the corresponding processors.
 */
int PMMG2D_distribute_mesh( PMMG2D_pParMesh parmesh )
{
  idx_t      *part;
  int        ier;

  // There is nothing to distribute on just 1 proc
  if( parmesh->nprocs == 1 ) return 1;

  // Partitions the mesh on root processor
  if( parmesh->myrank == parmesh->info.root ) {
    // Attribute to ref the point index
    for ( int i = 1; i <= parmesh->mesh->np; i++) parmesh->mesh->point[i].ref = i;

#ifdef USE_METIS
    // Call metis for partionning
    PMMG2D_CALLOC ( parmesh,part,parmesh->mesh->nt,idx_t,"allocate metis buffer", ier=5 );

    // Call metis, or recover a custom partitioning if provided (only to debug
    // the interface displacement, adaptation will be blocked) 
    if( !PMMG2D_PREDEF_PART ) {
      if ( PMMG2D_part_meshElts2metis( parmesh, part, parmesh->nprocs ) ) {
        // Store the partition number of each triangle in base
        for (int k = 1; k <= parmesh->mesh->nt; k++)
          parmesh->mesh->tria[k].base = part[k-1];
      }
      else ier = 5;
    } else {
      for(int k = 1; k <= parmesh->mesh->nt; k++ ) {
        part[k-1] = parmesh->mesh->tria[k].ref;
        parmesh->mesh->tria[k].tag[0] |= MG_REQ;
        parmesh->mesh->tria[k].tag[1] |= MG_REQ;
        parmesh->mesh->tria[k].tag[2] |= MG_REQ;
      }
    }

    PMMG2D_DEL_MEM(parmesh,part,idx_t,"deallocate metis buffer");
#else
    // Call scotch for partitioning
    ier = PMMG2D_Scotch_decomposition_root( parmesh );
#endif
  }

  // Send the local meshes from root process to the corresponding procs
  ier = PMMG2D_exchange_mesh( parmesh );

  int kk = MPI_Barrier(parmesh->comm);

  // Update the boundary entities in view of the future exchanges between processors
  // and to freeze the parallel fictitious boundaries during the remeshing process
  if ( !PMMG2D_mark_parallel_interface_nodes( parmesh ) ) {
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_LOWFAILURE);
  }

  // To debug, write one vtk file per processor containing the local mesh, metric and
  // parallel fictitious boundary points
  // check_vtk_file(parmesh);

  return ier;
}

