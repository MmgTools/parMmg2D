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
 * \file frontadvancing_pmmg2d.c
 * \brief Front advancing after a remeshing step
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 */

#include "parmmg2d.h"
#include "mpitypes_pmmg2d.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param list_index_pt number of triangles to exchange for each proc.
 * \param list_pt list of the triangles to exchange.
 * \param list_index_pp number of points to exchange potentially for each proc.
 * \param list_pp list of the points to exchange potentially.
 * \param list_mm list of the metric values to exchange potentially.
 * \param pt_global list of triangles received.
 * \param pp_global list of points received.
 * \param mm_global list of metric values received.
 * \param npt_global number of received triangles/
 * \param npp_global number of received points.
 * \return 0 if it fails, 1 otherwise.
 *
 * MPI exchange of the triangles between processors.
 *
 */
int PMMG2D_exchange( PMMG2D_pParMesh parmesh, int size, int* list_index_pt, MMG5_Tria** list_pt, 
                     int* list_index_pp, MMG5_Point** list_pp, double** list_mm,
                     MMG5_pTria* pt_global, MMG5_pPoint* pp_global, double** mm_global, int* npt_global, int* npp_global )
{

  int n, npt, npp;
  int ier = 1;
  int size_met = parmesh->met->size;
  MPI_Datatype   MPI_POINT, MPI_TRIANGLE;
  MPI_Request send_request, recv_request;
  MPI_Request send_request2, recv_request2;
  MPI_Status     status, status2;
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;

  *npt_global = 0;
  *npp_global = 0;

  PMMG2D_create_MPI_Tria( &MPI_TRIANGLE );
  PMMG2D_create_MPI_Point( &MPI_POINT );

  for ( n = 0; n < parmesh->nprocs; n++ ) {

    if (parmesh->myrank == n) continue;

    // Triangles
    MPI_CHECK( MPI_Send(&list_index_pt[n], 1, MPI_INT, n, (parmesh->myrank+1)*MPI_ANALYS_TAG + n, parmesh->comm), ier = 0 );

    MPI_CHECK( MPI_Recv(&npt, 1, MPI_INT, n, (n+1)*MPI_ANALYS_TAG + parmesh->myrank,
                         parmesh->comm, &status), ier = 0 );

    MPI_CHECK( MPI_Isend(list_pt[n], list_index_pt[n], MPI_TRIANGLE, n,
                         (parmesh->myrank+1)*MPI_ANALYS_TAG + parmesh->nprocs + 1 + n, parmesh->comm, &send_request), ier = 0 );

    *npt_global += npt;
    MMG5_pTria pt_local = (MMG5_pTria) malloc(npt * sizeof(MMG5_Tria));
    MPI_CHECK( MPI_Irecv(pt_local, npt, MPI_TRIANGLE, n,
                         (n+1)*MPI_ANALYS_TAG + parmesh->myrank + parmesh->nprocs + 1,
                         parmesh->comm, &recv_request), ier = 0 );
    MPI_Wait(&send_request, &status);
    MPI_Wait(&recv_request, &status);

    if (npt) {
      PMMG2D_REALLOC(mesh, *pt_global, *npt_global, *npt_global-npt,
                     MMG5_Tria,"Front advancing exchange triangle",return 0);
      memcpy(*pt_global + *npt_global - npt, pt_local, npt*sizeof(MMG5_Tria));
    }

    free(pt_local);

    // Points & Metrics
    MPI_CHECK( MPI_Send(&list_index_pp[n], 1, MPI_INT, n, (parmesh->myrank+1)*MPI_ANALYS_TAG + 2*parmesh->nprocs + 1 + n,
                         parmesh->comm), ier = 0 );

    MPI_CHECK( MPI_Recv(&npp, 1, MPI_INT, n, (n+1)*MPI_ANALYS_TAG + parmesh->myrank + 2*parmesh->nprocs + 1,
                         parmesh->comm, &status), ier = 0 );

    *npp_global += npp;

    MPI_CHECK( MPI_Isend(list_pp[n], list_index_pp[n], MPI_POINT, n,
                        (parmesh->myrank+1)*MPI_ANALYS_TAG + 3*parmesh->nprocs + 1 + n, parmesh->comm, &send_request), ier = 0 );

    MMG5_pPoint pp_local = (MMG5_pPoint) malloc(npp * sizeof(MMG5_Point));
    MPI_CHECK( MPI_Irecv(pp_local, npp, MPI_POINT, n,
                         (n+1)*MPI_ANALYS_TAG + parmesh->myrank + 3*parmesh->nprocs + 1,
                         parmesh->comm, &recv_request), ier = 0 );

    MPI_Wait(&send_request, &status);
    MPI_Wait(&recv_request, &status);

    MPI_CHECK( MPI_Isend(list_mm[n], size_met*list_index_pp[n], MPI_DOUBLE, n,
                        (parmesh->myrank+1)*MPI_ANALYS_TAG + 4*parmesh->nprocs + 1 + n, parmesh->comm, &send_request2), ier = 0 );

    MMG5_pSol mm_local = (MMG5_pSol) malloc(size_met*npp * sizeof(MMG5_Sol));
    MPI_CHECK( MPI_Irecv(mm_local, size_met*npp, MPI_DOUBLE, n,
                        (n+1)*MPI_ANALYS_TAG + parmesh->myrank + 4*parmesh->nprocs + 1,
                        parmesh->comm, &recv_request2), ier = 0 );

    MPI_Wait(&send_request2, &status2);
    MPI_Wait(&recv_request2, &status2);

    if (npp) {

      PMMG2D_REALLOC(mesh, *pp_global, *npp_global, *npp_global - npp,
                     MMG5_Point,"Front advancing exchange point",return 0);

      memcpy(*pp_global + *npp_global - npp, pp_local, npp*sizeof(MMG5_Point));

      PMMG2D_REALLOC(mesh, *mm_global, size_met* *npp_global, size_met* (*npp_global - npp),
                     double,"Front advancing exchange metric",return 0);

      memcpy(*mm_global + size_met*(*npp_global - npp), mm_local, size_met*npp*sizeof(double));

    }

    free(pp_local);
    free(mm_local);

  }

  PMMG2D_Free_MPI_meshDatatype(&MPI_POINT, &MPI_TRIANGLE, &MPI_TRIANGLE, false);

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if it fails, 1 otherwise.
 *
 * Remove the points that have been sent to another process.
 *
 */

int PMMG2D_remove_points( PMMG2D_pParMesh parmesh )
{
  int i,j,k;
  MMG5_pTria pt;
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;
  int* offset = (int*) malloc((mesh->np+1) * sizeof(int));
  int* list_points = (int*) malloc((mesh->np+1) * sizeof(int));

  for (i = 1; i <= mesh->np; i++) {
    list_points[i] = 0;
    offset[i] = 0;
  }

  // First store the list of points still in the mesh
  for ( k = 1; k <= mesh->nt; k++ ) {

    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for ( i = 0; i <= 2; i++ ) list_points[pt->v[i]] = 1;
  }

  // Second remove the points and update the metrics
  j = 0;
  for ( k = 1; k <= mesh->np; k++ ) {

    if (!list_points[k]) {
      j += 1;
      offset[k] = j;
      continue;
    }

    mesh->point[k-j] = mesh->point[k];
    offset[k] = j;

    if (parmesh->met->size == 1) { // isotropic
      parmesh->met->m[k-j] = parmesh->met->m[k];
    }
    else {
      for ( i = 0; i < 3; i++) { // anisotropic
        parmesh->met->m[3*(k-j)+i] = parmesh->met->m[3*k+i];
      }
    }

  }

  if ( !j ) {
    free(offset);
    free(list_points);
    return 1;
  }

  mesh->np -= j;

  PMMG2D_REALLOC(mesh, mesh->point, mesh->np+1, mesh->np+j+1,
                 MMG5_Point,"Front advancing realloc points", return 0);
  PMMG2D_REALLOC(mesh, parmesh->met->m, parmesh->met->size*(mesh->np+1), parmesh->met->size*(mesh->np+j+1),
                 double, "Front advancing realloc metrics", return 0);

  // Third correct the index number in the triangles
  for ( k = 1; k <= mesh->nt; k ++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    for ( i = 0; i <= 2; i++ ) mesh->tria[k].v[i] -= offset[mesh->tria[k].v[i]];

  }

  free(offset);
  free(list_points);

  return 1;

}

static int compare(const void *a, const void *b) {
  return (*(int*)a - *(int*)b);
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param list_index array containing the indices of the exchanged triangles.
 * \param list_pt list of the triangles to exchange.
 * \param list_index_pt number of triangles to exchange for each proc.
 * \param list_pp list of the points to exchange potentially.
 * \param list_index_pp number of points to exchange potentially for each proc.
 * \param list_mm list of the metric values to exchange potentially.
 * \param list_deleted_triangles list of deleted triangles.
 * \param size_list_deleted_triangles size of the list of deleted triangles.
 *
 * Add adjacent triangles in the list of exchanged interface triangles.
 *
 */
void PMMG2D_add_layer_triangles( PMMG2D_pParMesh parmesh, int** list_index, 
                                 MMG5_Tria** list_pt,
				 int* list_index_pt, 
                                 MMG5_Point** list_pp,
				 int* list_index_pp, 
                                 double** list_mm,
				 int* list_deleted_triangles,
                                 int size_list_deleted_triangles )
{
  int i, j, k, l, m, n, p, k1;
  MMG5_pTria pt;
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;
  int start[parmesh->nprocs];

  p = size_list_deleted_triangles;

  if (p == 0) return; // No triangle to send

  for (i = 0; i < parmesh->nprocs; i++) start[i] = 0;

  int** list_extra_triangles = (int**) malloc(parmesh->nprocs * sizeof(int*));

  for (i = 1; i < parmesh->info.ifc_layers; i++) { // Loop over the layers of triangles

    for (j = 0; j < parmesh->nprocs; j++) list_extra_triangles[j] = (int*) malloc(sizeof(int));
    int size[parmesh->nprocs];

    for (j = 0; j < parmesh->nprocs; j++) size[j] = 0;

    for (j = 0; j < parmesh->nprocs; j++) { // Loop over the processors

      if ( parmesh->myrank == j ) continue;

      for (l = start[j]; l < list_index_pt[j]; l++) { // Loop over the triangles to exchange with this proc

        k = list_index[j][l];

        for( n = 1; n < 4; n++ ) { // Loop over the adjacent triangles to triangle k
  
          k1 = mesh->adja[3*(k-1)+n];
          if( !k1 ) continue;
          k1 /= 3;
          pt = &mesh->tria[k1];

          if ( pt->flag ) continue; // The triangle is already exchanged

          list_extra_triangles[j][size[j]++] = k1;
          pt->flag = mesh->base;
          list_extra_triangles[j] = (int*) realloc(list_extra_triangles[j], (size[j]+1) * sizeof(int));
          
        }

      }
    
    }

    // Add the extra triangles to the lists
    for (j = 0; j < parmesh->nprocs; j++) {

      start[j] = list_index_pt[j];

      for (l = 0; l < size[j]; l++) {
    
        k = list_extra_triangles[j][l];

        pt = &mesh->tria[k];

        list_index[j][list_index_pt[j]] = k;

        list_deleted_triangles[p] = k;

        p++;

        list_pt[j][list_index_pt[j]] = mesh->tria[k];
        list_index_pt[j]++;

	list_pt[j] = (MMG5_pTria) realloc(list_pt[j], (list_index_pt[j]+1) * sizeof(MMG5_Tria));
        list_index[j] = (int*) realloc(list_index[j], (list_index_pt[j]+1) * sizeof(int));

        for (n = 0; n <= 2; n++) {
          list_pp[j][list_index_pp[j]] = mesh->point[pt->v[n]];
          if (parmesh->met->size == 1) { // isotropic
            list_mm[j][list_index_pp[j]] = parmesh->met->m[pt->v[n]];
          }
          else { // anisotropic
            list_mm[j][3*list_index_pp[j]] = parmesh->met->m[3*pt->v[n]];
            list_mm[j][3*list_index_pp[j]+1] = parmesh->met->m[3*pt->v[n]+1];
            list_mm[j][3*list_index_pp[j]+2] = parmesh->met->m[3*pt->v[n]+2];
          }
          list_index_pp[j]++;
        }

        list_pp[j] = (MMG5_pPoint) realloc(list_pp[j], (list_index_pp[j]+3) * sizeof(MMG5_Point));
        list_mm[j] = (double*) realloc(list_mm[j], parmesh->met->size*(list_index_pp[j]+3) * sizeof(double));
      }

    }

    for (j = 0; j < parmesh->nprocs; j++) free(list_extra_triangles[j]);

  }

  free(list_extra_triangles);

  // Sort list_deleted_triangles
  qsort(list_deleted_triangles, p, sizeof(int), compare);

}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param list_index_pt number of triangles to exchange for each proc.
 * \param list_pt list of the triangles to exchange.
 * \param list_index_pp number of points to exchange potentially for each proc.
 * \param list_pp list of the points to exchange potentially.
 * \param list_mm list of the metric values to exchange potentially.
 * \param list_deleted_triangles list of deleted triangles.
 * \return 0 if it fails, 1 otherwise.
 *
 * Find and store the interface triangles to exchange between processors.
 *
 */
int PMMG2D_find_exchanged_triangles( PMMG2D_pParMesh parmesh, int* list_index_pt, MMG5_Tria** list_pt,
                                     int* list_index_pp, MMG5_Point** list_pp,
                                     double** list_mm,
				     int* list_deleted_triangles )
{
  int nb_triangles_proc[parmesh->nprocs];
  MMG5_pMesh mesh;
  MMG5_pTria pt;
  mesh = parmesh->mesh;
  int is_interface_tria, min_nt, send_proc, i, j, k, l, n;
  int* inverse_interface_node_list = (int*) malloc((mesh->np+1) * sizeof(int));
  int** list_index = (int**) malloc(parmesh->nprocs * sizeof(int*));
  for (i = 0; i < parmesh->nprocs; i++) list_index[i] = (int*) malloc(sizeof(int));

  // Initialisation
  for (i = 0; i < parmesh->nprocs; i++) {
    list_index_pt[i] = 0;
    list_index_pp[i] = 0;
  }
  for (i = 0; i < mesh->nt; i++) list_deleted_triangles[i] = 0;
  for (i = 0; i <= mesh->np; i++) inverse_interface_node_list[i] = 0;
  for (i = 0; i < parmesh->nip; i++) 
    inverse_interface_node_list[parmesh->list_interface_nodes[i].index[0]] = i;

  for ( k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    pt->flag = 0;
  }

  mesh->base++;

  // Exchange the number of triangles per process
  MPI_CHECK( MPI_Allgather(&mesh->nt, 1, MPI_INT, nb_triangles_proc, 1, MPI_INT, parmesh->comm),
             return 0 );

  // Find the triangles to exchange
  n = 0;
  for (k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    // Check if the triangle gets at least one vertex on the interface between two processors
    is_interface_tria = 0;
    for (i = 0; i <= 2; i++) {
      if (mesh->point[pt->v[i]].tag & MG_PARBDY) {
        is_interface_tria = 1;
        break;
      }
    }

    if (!is_interface_tria) continue;

    // Find the adjacent processors that contain the smallest number of triangles
    min_nt = mesh->nt;
    send_proc = parmesh->myrank;

    for (i = 0; i <= 2; i++) {
      if (mesh->point[pt->v[i]].tag & MG_PARBDY) {
        l = inverse_interface_node_list[pt->v[i]];
        for (j = 1; j < parmesh->list_interface_nodes[l].nb; j++) {
          if (nb_triangles_proc[parmesh->list_interface_nodes[l].proc[j]] < min_nt) {
            min_nt = nb_triangles_proc[parmesh->list_interface_nodes[l].proc[j]];
            send_proc = parmesh->list_interface_nodes[l].proc[j];
          }
        }
      }
    }

    // Check if the triangle must remain in the own process or if another neighboring
    // processor contains less triangles
    if (send_proc == parmesh->myrank) continue;

    list_deleted_triangles[n] = k;
    pt->flag = mesh->base;

    n++;

    list_pt[send_proc][list_index_pt[send_proc]] = mesh->tria[k];
    list_index[send_proc][list_index_pt[send_proc]] = k;
    list_index_pt[send_proc]++;

    list_pt[send_proc] = (MMG5_pTria) realloc(list_pt[send_proc], (list_index_pt[send_proc]+1) * sizeof(MMG5_Tria));
    list_index[send_proc] = (int*) realloc(list_index[send_proc], (list_index_pt[send_proc]+1) * sizeof(int));

    for (i = 0; i <= 2; i++) {
      list_pp[send_proc][list_index_pp[send_proc]] = mesh->point[pt->v[i]];
      if (parmesh->met->size == 1) { // isotropic
        list_mm[send_proc][list_index_pp[send_proc]] = parmesh->met->m[pt->v[i]];
      }
      else { // anisotropic
        list_mm[send_proc][3*list_index_pp[send_proc]] = parmesh->met->m[3*pt->v[i]];
        list_mm[send_proc][3*list_index_pp[send_proc]+1] = parmesh->met->m[3*pt->v[i]+1];
        list_mm[send_proc][3*list_index_pp[send_proc]+2] = parmesh->met->m[3*pt->v[i]+2];
      }
      list_index_pp[send_proc]++;
    }

    list_pp[send_proc] = (MMG5_pPoint) realloc(list_pp[send_proc], (list_index_pp[send_proc]+3) * sizeof(MMG5_Point));
    list_mm[send_proc] = (double*) realloc(list_mm[send_proc], parmesh->met->size*(list_index_pp[send_proc]+3) * sizeof(double));
  }

  free(inverse_interface_node_list);

  // Potentially add a layer of triangles in addition to the interface triangles
  PMMG2D_add_layer_triangles( parmesh, list_index, list_pt, list_index_pt, list_pp, 
                              list_index_pp, list_mm, list_deleted_triangles, n );

  for (i = 0; i < parmesh->nprocs; i++) free(list_index[i]);
  free(list_index);

  return 1;
  
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param pt_global list of triangles received.
 * \param pp_global list of points received.
 * \param mm_global list of metric values received.
 * \param npt_global number of received triangles/
 * \param npp_global number of received points.
 * \return 0 if it fails, 1 otherwise.
 *
 * Add the received triangles to the local meshes.
 *
 */
int PMMG2D_add_triangles( PMMG2D_pParMesh parmesh, MMG5_pTria pt_global, 
                          MMG5_pPoint pp_global, double* mm_global, int npt, int npp )
{
  int i, j, k, max_glob;
  MMG5_pMesh mesh;
  MMG5_pSol met;
  mesh = parmesh->mesh;
  met = parmesh->met;

  int max = 0;
  for (k = 1; k <= mesh->np; k++) {
    if (mesh->point[k].ref > max) max = mesh->point[k].ref;
  }

  for (k = 0; k < npp; k++) {
    if (pp_global[k].ref > max) max = pp_global[k].ref;
  }

  MPI_Allreduce( &max, &max_glob, 1, MPI_INT, MPI_MAX, parmesh->comm );
  int* list_ref = (int*) malloc((max_glob+1)*sizeof(int));

  for (k = 0; k <= max_glob; k++) list_ref[k] = 0;
  for (k = 1; k <= mesh->np; k++) list_ref[mesh->point[k].ref] = k;

  j = 0;
  for (k = 0; k < npt; k++) {
    for (i = 0; i <= 2; i++) {

      // Check whether the point already exists in the local mesh
      // We use the ref number which is unique all across the global mesh
      if (list_ref[pp_global[j].ref] && j < npp) {
        if (mesh->point[list_ref[pp_global[j].ref]].tag < pp_global[j].tag ) mesh->point[list_ref[pp_global[j].ref]].tag = pp_global[j].tag;
        if (mesh->point[list_ref[pp_global[j].ref]].xp < pp_global[j].xp ) mesh->point[list_ref[pp_global[j].ref]].xp = pp_global[j].xp;
        pt_global[k].v[i] = list_ref[pp_global[j].ref];
      }
      else {
        mesh->np++;
        met->np++;
        PMMG2D_REALLOC(mesh, mesh->point, mesh->np+1, mesh->np,
                       MMG5_Point,"Add one point in front advancement", return 0);
        memcpy(mesh->point + mesh->np, &pp_global[j], sizeof(MMG5_Point));
        list_ref[pp_global[j].ref] = mesh->np;
        pt_global[k].v[i] = mesh->np;

        PMMG2D_REALLOC(mesh, met->m, met->size*(mesh->np+1), met->size*mesh->np,
                       double,"Add one point in front advancement", return 0);
        memcpy(met->m + met->size*mesh->np, &mm_global[met->size*j], met->size*sizeof(double));

      }
      j += 1;
    }
  }

  free(list_ref);

  PMMG2D_REALLOC(mesh, mesh->tria, mesh->nt+1+npt, mesh->nt+1,
                 MMG5_Tria,"Front advancing realloc triangles", return 0);
  memcpy(mesh->tria + mesh->nt + 1, pt_global, npt*sizeof(MMG5_Tria));
  mesh->nt += npt;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param list_adjacent_triangles lists of adjacent triangles.
 * \param number_adjacent_triangles number of adjacent triangles in the lists.
 *
 * Build lists of adjacent triangles. If there is no isolated triangle, there is only one list.
 *
 */
void PMMG2D_adjacent_triangles( MMG5_pMesh mesh, int *list_adjacent_triangles, int* number_adjacent_triangles )
{
  int i,j,k,k1,l,base,size_list,start,cur;
  MMG5_pTria pt, pt1;

  for ( k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    pt->flag = 0;
  }

  // Build lists of adjacent triangles
  base = ++mesh->base;
  i = 0;
  j = 1;
  cur = 0;

  while (j <= mesh->nt) {

    // Find a non-visited starting triangle
    start = -1;
    for ( k = 1; k <= mesh->nt; k++) {
      pt = &mesh->tria[k];
      if ( !MG_EOK(pt) ) continue;
      if ( !pt->flag ) {
        start = k;
        break;
      }
    }

    if (start == -1) fprintf(stderr, "ERROR in PMMG2D_adjacent_triangles\n");

    // Store initial triangle
    list_adjacent_triangles[cur] = start;

    size_list = 1;

    // Flag initial triangle
    mesh->tria[start].flag = base;

    // Explore list and fill it by adjacency
    while( cur < size_list + j - 1 ) {
      k = list_adjacent_triangles[cur];

      // Loop on adjacent triangles
      for( l = 1; l < 4; l++ ) {
        k1 = mesh->adja[3*(k-1)+l];
        if( !k1 ) continue;

        k1 /= 3;
        pt1 = &mesh->tria[k1];

        // Skip already visited triangles (by this or another list)
        if( !pt1->flag ) {
          list_adjacent_triangles[j - 1 + size_list++] = k1;
          pt1->flag = base;
        }
      }

      cur++;
    }

    j += size_list;

    number_adjacent_triangles[i] = size_list;

    i += 1;

  }

}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param list_adjacent_triangles lists of adjacent triangles.
 * \param number_adjacent_triangles number of adjacent triangles in the lists.
 * \param proc_adjacent_triangles processor number adjacent to at least one triangle of the list
 * \param n_list number of distinct list of adjacent triangles in the mesh
 * \param n_longest index number of the longest list
 *
 * Build lists of adjacent triangles. If there is no isolated triangle, there is only one list.
 *
 */
void PMMG2D_find_adjacent_processor( PMMG2D_pParMesh parmesh, int** list_adjacent_triangles, 
                                     int* number_adjacent_triangles, int* proc_adjacent_triangles, 
                                     int n_list, int n_longest )
{
  int i, j, l, max_ref, must_continue;
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;
  MMG5_pTria pt;

  int* list_ref = (int*) malloc(mesh->np * sizeof(int));
  int rcounts[parmesh->nprocs];
  int displs[parmesh->nprocs];

  // List the global references on each proc
  max_ref = 0;
  for (i = 1; i <= mesh->np; i++) {
    list_ref[i-1] = mesh->point[i].ref;
    if (mesh->point[i].ref > max_ref) max_ref = mesh->point[i].ref;
  }

  MPI_Allgather(&mesh->np, 1, MPI_INT, rcounts, 1, MPI_INT, parmesh->comm);

  displs[0] = 0;
  for (i = 1; i < parmesh->nprocs; i++) displs[i] = displs[i-1] + rcounts[i-1];

  // Gather the list of global references
  int* list_global_ref = (int*) malloc((displs[parmesh->nprocs-1]+rcounts[parmesh->nprocs-1]) * sizeof(int));
  MPI_Allgatherv(&list_ref[0], mesh->np, MPI_INT, list_global_ref, rcounts, displs, MPI_INT, parmesh->comm);

  free(list_ref);

  int max_ref_global;
  MPI_Allreduce( &max_ref, &max_ref_global, 1, MPI_INT, MPI_MAX, parmesh->comm );

  int* global_ref_proc = (int*) malloc((max_ref_global+1) * sizeof(int));
  for (i = 0; i <= max_ref_global; i++) global_ref_proc[i] = -1;

  int num_proc = 0;
  for (i = 0;  i < displs[parmesh->nprocs-1]+rcounts[parmesh->nprocs-1]; i++) {
    if ( num_proc < parmesh->nprocs-1 && i == displs[num_proc+1] ) {
      num_proc++;
      while ( num_proc < parmesh->nprocs-1 && rcounts[num_proc] == 0 ) num_proc++;
    }
    if ( num_proc == parmesh->myrank ) continue;
    global_ref_proc[list_global_ref[i]] = num_proc;
  }

  // Attribute to each list of adjacent triangles an adjacent processor
  for (i = 0; i < n_list; i++) {

    if (i == n_longest) continue;
    
    must_continue = 1;

    for (j = 0; j < number_adjacent_triangles[i]; j++) {
      pt = &mesh->tria[list_adjacent_triangles[i][j]];
      if ( !MG_EOK(pt) ) continue;

      for (l = 0; l <= 2; l++) {
        if ( global_ref_proc[mesh->point[pt->v[l]].ref] == -1 ) continue;
        proc_adjacent_triangles[i] = global_ref_proc[mesh->point[pt->v[l]].ref];
        must_continue = 0;
        break;
      }

      if ( !must_continue ) break;

    }

  }

  free(list_global_ref);
  free(global_ref_proc);

}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param nt_tmp integer to create arrays with a smaller size than the number of points.
 * \return 0 if it fails, 1 otherwise. 
 *
 * Check the contiguity of the local mesh and exchange triangles 
 * to avoid the existence of isolated triangles. 
 *
 */
int PMMG2D_check_contiguity( PMMG2D_pParMesh parmesh, int nt_tmp )
{
  MMG5_pTria pt;
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;
  int* list_deleted_triangles = (int*) malloc((mesh->nt+1) * sizeof(int));
  int* list_adjacent_triangles = (int*) malloc(mesh->nt * sizeof(int));
  int* number_adjacent_triangles = (int*) malloc((mesh->nt - nt_tmp + 1) * sizeof(int));
  int list_index_pt[parmesh->nprocs];
  int list_index_pp[parmesh->nprocs];
  MMG5_pTria pt_global = NULL;
  MMG5_pPoint pp_global = NULL;
  double* mm_global = NULL;
  int i, j, k, n_list, n_longest, max, send_proc, npt_global, npp_global, ier, ieresult;

  // 1) Compute the mesh adjacency
  PMMG2D_DEL_MEM(mesh, mesh->adja, int, "remove adjacency table" );

  if ( (!mesh->adja) && (mesh->nt) && (MMG2D_hashTria( mesh ) != 1) ) {
    fprintf( stderr,"  ## PMMG2D Hashing problem (2).\n" );
    return 0;
  }

  assert(mesh->nt >= nt_tmp);

  // 2) Build lists of adjacent triangles
  for ( i = 0; i <= mesh->nt - nt_tmp; i++ ) number_adjacent_triangles[i] = 0;
  PMMG2D_adjacent_triangles( mesh, list_adjacent_triangles, number_adjacent_triangles );

  // Compute the number of lists and find the longest one
  n_list = 1;
  n_longest = 0;
  max = number_adjacent_triangles[0];
  while ( n_list <= mesh->nt - nt_tmp && number_adjacent_triangles[n_list] ) {
    if ( number_adjacent_triangles[n_list] > max ) {
      n_longest = n_list;
      max = number_adjacent_triangles[n_list];
    }
    n_list++;
  }

  int** list_adjacent_triangles_decomposed = (int**) malloc(n_list * sizeof(int*));
  for (i = 0; i < n_list; i++) list_adjacent_triangles_decomposed[i] = (int*) malloc(max * sizeof(int));
  k = 0;
  for ( i = 0; i < n_list; i++ ) {
    for ( j = 0; j < number_adjacent_triangles[i]; j++ ) {
      list_adjacent_triangles_decomposed[i][j] = list_adjacent_triangles[k];
      k++;
    }
  }

  free(list_adjacent_triangles);

  // 3) Find an adjacent processor for each list except the longest one
  int proc_adjacent_triangles[n_list];
  PMMG2D_find_adjacent_processor( parmesh, list_adjacent_triangles_decomposed, number_adjacent_triangles, 
                                  proc_adjacent_triangles, n_list, n_longest );

  // 4) Loop over the lists (except the longest one) and 
  //    create the lists of points, triangles and metrics to exchange
  for ( i = 0; i < parmesh->nprocs; i++ ) {
    list_index_pt[i] = 0;
    list_index_pp[i] = 0;
  }
  for (i = 0; i <= mesh->nt; i++) list_deleted_triangles[i] = 0;

  int sum_adjacent_triangles = 0;
  for (j = 0; j < n_list; j++) {
    if ( j == n_longest ) continue;
    sum_adjacent_triangles += number_adjacent_triangles[j];
  }

  MMG5_pTria* list_pt = (MMG5_pTria*) malloc(parmesh->nprocs * sizeof(MMG5_pTria));
  for (k = 0; k < parmesh->nprocs; k++) list_pt[k] = (MMG5_pTria) malloc(sum_adjacent_triangles * sizeof(MMG5_Tria));
  MMG5_pPoint* list_pp = (MMG5_pPoint*) malloc(parmesh->nprocs * sizeof(MMG5_pPoint));
  for (k = 0; k < parmesh->nprocs; k++) list_pp[k] = (MMG5_pPoint) malloc(3*sum_adjacent_triangles * sizeof(MMG5_Point));
  double** list_mm = (double**) malloc(parmesh->nprocs * sizeof(double*));
  for (k = 0; k < parmesh->nprocs; k++) list_mm[k] = (double*) malloc(parmesh->met->size*3*sum_adjacent_triangles * sizeof(double));

  for ( j = 0; j < n_list; j++) {

    if ( j == n_longest ) continue;

    send_proc = proc_adjacent_triangles[j];

    for ( k = 0; k < number_adjacent_triangles[j]; k++ ) {

      pt = &mesh->tria[list_adjacent_triangles_decomposed[j][k]];
      list_deleted_triangles[list_adjacent_triangles_decomposed[j][k]] = 1; 

      list_pt[send_proc][list_index_pt[send_proc]] = mesh->tria[list_adjacent_triangles_decomposed[j][k]];

      list_index_pt[send_proc]++;

      for (i = 0; i <= 2; i++) {

        list_pp[send_proc][list_index_pp[send_proc]] = mesh->point[pt->v[i]];

        if (parmesh->met->size == 1) { // isotropic
          list_mm[send_proc][list_index_pp[send_proc]] = parmesh->met->m[pt->v[i]];
        }
        else { // anisotropic
          list_mm[send_proc][3*list_index_pp[send_proc]] = parmesh->met->m[3*pt->v[i]];
          list_mm[send_proc][3*list_index_pp[send_proc]+1] = parmesh->met->m[3*pt->v[i]+1];
          list_mm[send_proc][3*list_index_pp[send_proc]+2] = parmesh->met->m[3*pt->v[i]+2];
        }

        list_index_pp[send_proc]++;

      }
    }

  }

  free(number_adjacent_triangles);
  for (i = 0; i < n_list; i++) free(list_adjacent_triangles_decomposed[i]);
  free(list_adjacent_triangles_decomposed);

  // 5) Exchange the triangles, points and metrics
  ier = PMMG2D_exchange( parmesh, sum_adjacent_triangles, list_index_pt, list_pt, list_index_pp, list_pp, list_mm,
                         &pt_global, &pp_global, &mm_global, &npt_global, &npp_global );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error in the mpi exchanges during the contiguity check. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  for (k = 0; k < parmesh->nprocs; k++) {
    free(list_pt[k]);
    free(list_pp[k]);
    free(list_mm[k]);
  }
  free(list_pt);
  free(list_pp);
  free(list_mm);

  // 6) Remove the sent triangles
  j = 0;
  for ( k = 1; k <= mesh->nt; k++ ) {
    if ( list_deleted_triangles[k] ) {
      j += 1;
      continue;
    }
    mesh->tria[k-j] = mesh->tria[k];
  }

  free(list_deleted_triangles);

  mesh->nt -= j;
  PMMG2D_REALLOC(mesh, mesh->tria, mesh->nt+1, mesh->nt+j+1,
                 MMG5_Tria,"contiguity check realloc triangles",return 0);

  // 7) Remove the points that do not belong to the mesh anymore
  ier = PMMG2D_remove_points( parmesh );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when removing the points during the contiguity check. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 8) Add the triangles to their corresponding new process
  ier = PMMG2D_add_triangles(parmesh, pt_global, pp_global, mm_global, npt_global, npp_global);
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when adding the new triangles during the contiguity check. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Free the pointers
  PMMG2D_DEL_MEM(mesh, pt_global, MMG5_Tria, "Delete pt_global in front advancing process");
  PMMG2D_DEL_MEM(mesh, pp_global, MMG5_Point, "Delete pp_global in front advancing process");
  PMMG2D_DEL_MEM(mesh, mm_global, double, "Delete mm_global in front advancing process");

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * Untag parallel nodes and edges.
 *
 */
void PMMG2D_untag_parallel( PMMG2D_pParMesh parmesh )
{
  int i,k;
  MMG5_pTria pt;
  MMG5_pPoint ppt;

  // Parallel edges
  for (k = 1; k <= parmesh->mesh->nt; k++) {
    pt = &parmesh->mesh->tria[k];
  
    if ( !MG_EOK(pt) ) continue;
  
    for (i = 0; i <= 2; i++) {
      if ( pt->tag[i] & MG_PARBDY ) {
        pt->tag[i] &= ~MG_PARBDY;
        if (!(pt->tag[i] & MG_PARBDYBDY)) pt->tag[i] &= ~MG_REQ;
        if (!(pt->tag[i] & MG_REQ) || !(pt->tag[i] & MG_BDY)) pt->tag[i] &= ~MG_NOSURF;
        pt->tag[i] &= ~MG_PARBDYBDY;
      }
    }
  }

  // Parallel points
  for (k = 1; k <= parmesh->mesh->np; k++) {
    ppt = &parmesh->mesh->point[k];

    if ( ppt->tag & MG_PARBDY ) {
      ppt->tag &= ~MG_PARBDY;
      if (!(ppt->tag & MG_PARBDYBDY)) ppt->tag &= ~MG_REQ;
      if (!(ppt->tag & MG_REQ) || !(ppt->tag & MG_BDY)) ppt->tag &= ~MG_NOSURF;
      ppt->tag &= ~MG_PARBDYBDY;
    }
  }

}

void PMMG2D_update( PMMG2D_pParMesh parmesh )
{
  parmesh->mesh->npi = parmesh->mesh->np;
  parmesh->mesh->npnil = 0;
  parmesh->met->np = parmesh->mesh->np;
  parmesh->met->npi = parmesh->mesh->np;
  parmesh->mesh->npmax = parmesh->mesh->np;
  parmesh->met->npmax = parmesh->met->np;

  parmesh->mesh->nenil = 0;
  parmesh->mesh->ntmax = parmesh->mesh->nti = parmesh->mesh->nt;
  parmesh->mesh->xtmax = parmesh->mesh->xt;

}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if it fails, 1 otherwise. 
 *
 * Modify the interface boundaries between processors
 *
 */
int PMMG2D_frontadvancing( PMMG2D_pParMesh parmesh )
{
  int j, k, npt_global, npp_global, ier, ieresult, nt_tmp;
  int list_index_pt[parmesh->nprocs];
  int list_index_pp[parmesh->nprocs];
  MMG5_pMesh mesh;
  MMG5_pTria pt_global = NULL;
  MMG5_pPoint pp_global = NULL;
  double* mm_global = NULL;
  mesh = parmesh->mesh;
  MMG5_pTria* list_pt = (MMG5_pTria*) malloc(parmesh->nprocs * sizeof(MMG5_pTria)); 
  for (k = 0; k < parmesh->nprocs; k++) list_pt[k] = (MMG5_pTria) malloc(sizeof(MMG5_Tria));
  MMG5_pPoint* list_pp = (MMG5_pPoint*) malloc(parmesh->nprocs * sizeof(MMG5_pPoint)); 
  for (k = 0; k < parmesh->nprocs; k++) list_pp[k] = (MMG5_pPoint) malloc(3*sizeof(MMG5_Point));
  double** list_mm = (double**) malloc(parmesh->nprocs * sizeof(double*)); 
  for (k = 0; k < parmesh->nprocs; k++) list_mm[k] = (double*) malloc(3*parmesh->met->size*sizeof(double));
  int* list_deleted_triangles = (int*) malloc((mesh->nt+1)*sizeof(int));

  // 1) Find the triangles to exchange
  ier = PMMG2D_find_exchanged_triangles( parmesh, list_index_pt, list_pt, list_index_pp, list_pp, list_mm, list_deleted_triangles );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when searching the triangles to exchange during the front advancement. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 2) Exchange the triangles, points and metrics
  ier = PMMG2D_exchange( parmesh, mesh->nt, list_index_pt, list_pt, list_index_pp, list_pp, list_mm, 
                         &pt_global, &pp_global, &mm_global, &npt_global, &npp_global );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error in the mpi exchanges during the front advancement. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  for (k = 0; k < parmesh->nprocs; k++) {
    free(list_pt[k]);
    free(list_pp[k]);
    free(list_mm[k]);
  }
  free(list_pt);
  free(list_pp);
  free(list_mm);

  // 3) Remove the sent triangles from the local mesh
  j = 0;
  for ( k = 1; k <= mesh->nt; k++ ) {
    if (list_deleted_triangles[j] == k) {
      j += 1;
      continue;
    }
    mesh->tria[k-j] = mesh->tria[k];
  }

  mesh->nt -= j;
  nt_tmp = mesh->nt - j; // temporary number of triangles, used to check the contiguity (removing twice the triangles)

  free(list_deleted_triangles);

  PMMG2D_REALLOC(mesh, mesh->tria, mesh->nt+1, mesh->ntmax+1,
                 MMG5_Tria,"Front advancing realloc triangles",return 0);
  PMMG2D_REALLOC(mesh, mesh->point, mesh->np+1, mesh->npmax+1,
                 MMG5_Point,"Front advancing realloc points", return 0);
  PMMG2D_REALLOC(mesh, parmesh->met->m, parmesh->met->size*(mesh->np+1), parmesh->met->size*(mesh->npmax+1),
                 double, "Front advancing realloc metrics", return 0);

  // 4) Remove the points that do not belong to the mesh anymore
  ier = PMMG2D_remove_points( parmesh ); 
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when removing the points during the front advancement. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 5) Add the triangles to their corresponding new process
  ier = PMMG2D_add_triangles( parmesh, pt_global, pp_global, mm_global, npt_global, npp_global);
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when adding the new triangles during the front advancement. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  PMMG2D_update( parmesh );

  // 6) Check the contiguity and correct the isolated triangles
  ier = PMMG2D_check_contiguity( parmesh, nt_tmp );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when checking the contiguity during the front advancement. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  PMMG2D_update( parmesh );

  // 7) Untag the interface parallel points and edges
  PMMG2D_untag_parallel( parmesh );

  // 8) Mark the interface nodes
  ier = PMMG2D_mark_parallel_interface_nodes( parmesh );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when marking the interface nodes during the contiguity check. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 9) Update the list of interface nodes
  PMMG2D_free_interface_nodes_list( parmesh );
  PMMG2D_fill_interface_nodes_list( parmesh );

  // Free the pointers
  PMMG2D_DEL_MEM(mesh, pt_global, MMG5_Tria, "Delete pt_global in front advancing process");
  PMMG2D_DEL_MEM(mesh, pp_global, MMG5_Point, "Delete pp_global in front advancing process");
  PMMG2D_DEL_MEM(mesh, mm_global, double, "Delete mm_global in front advancing process");

  return 1;
}
