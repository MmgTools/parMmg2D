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
 * \file mergemesh_pmmg2d.c
 * \brief Merge the mesh.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 */
#include "parmmg2d.h"
#include "mpitypes_pmmg2d.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if fail, 1 if successc(on all procs)
 *
 *  merge all meshes and metrics to single mesh and metric in P0's parmesh
 */
int PMMG2D_merge_parmesh( PMMG2D_pParMesh parmesh ) 
{
  int n, k, i, i1, i2, ier;
  MMG5_pMesh mesh;
  MMG5_pSol met;
  MMG5_pTria pt;
  MPI_Datatype   MPI_POINT, MPI_TRIANGLE; 
  MPI_Status     status;
  int list_nt[parmesh->nprocs], list_np[parmesh->nprocs];
  int nt_glob = 0, np_glob = 0, nm_glob = 0;
  int offset = 0;

  mesh = parmesh->mesh;
  met = parmesh->met;

  // 1) Gather the mesh and metric information in root process
  for (n = 0; n < parmesh->nprocs; n++) {
    list_nt[n] = 0;
    list_np[n] = 0;
  }

  // Create specific types to MPI exchanges
  PMMG2D_create_MPI_Tria( &MPI_TRIANGLE );
  PMMG2D_create_MPI_Point( &MPI_POINT );

  // Triangles
  MPI_CHECK( MPI_Gather(&mesh->nt, 1, MPI_INT, list_nt, 1, MPI_INT, parmesh->info.root, parmesh->comm), ier = 0 );

  for ( n = 0; n < parmesh->nprocs; n++ ) nt_glob += list_nt[n];
  int displs[parmesh->nprocs];
  displs[0] = 0;
  for ( n = 1; n < parmesh->nprocs; n++ ) displs[n] = displs[n-1] + list_nt[n-1];

  // Problem when using MPI_Gather(v) with great meshes (probably due to MPI memory issues).
  // Send and Recv are therefore used so that MPI deals with smaller arrays each exchange.
  MMG5_pTria list_triangles = (MMG5_pTria) malloc(nt_glob * sizeof(MMG5_Tria));
  if (parmesh->myrank == parmesh->info.root) {
    for ( n = 0; n < mesh->nt; n++) list_triangles[n] = mesh->tria[n+1];
    for ( n = 1; n < parmesh->nprocs; n++) {
      MPI_CHECK( MPI_Recv( &list_triangles[displs[n]], list_nt[n], MPI_TRIANGLE, n, (n+1)*MPI_ANALYS_TAG, parmesh->comm, &status), ier = 0 );
    }
  }
  else MPI_CHECK( MPI_Send(&mesh->tria[1], mesh->nt, MPI_TRIANGLE, parmesh->info.root, (parmesh->myrank+1)*MPI_ANALYS_TAG, parmesh->comm), ier = 0 );

  // Points
  MPI_CHECK( MPI_Gather(&mesh->np, 1, MPI_INT, list_np, 1, MPI_INT, parmesh->info.root, parmesh->comm), ier = 0 );

  for ( n = 0; n < parmesh->nprocs; n++ ) np_glob += list_np[n];
  int displsp[parmesh->nprocs];
  displsp[0] = 0;
  for ( n = 1; n < parmesh->nprocs; n++ ) displsp[n] = displsp[n-1] + list_np[n-1];

  MMG5_pPoint list_points = (MMG5_pPoint) malloc(np_glob * sizeof(MMG5_Point));
  if (parmesh->myrank == parmesh->info.root) {
    for ( n = 0; n < mesh->np; n++) list_points[n] = mesh->point[n+1];
    for ( n = 1; n < parmesh->nprocs; n++) {
      MPI_CHECK( MPI_Recv( &list_points[displsp[n]], list_np[n], MPI_POINT, n, (n+1)*MPI_ANALYS_TAG + parmesh->nprocs + 1, parmesh->comm, &status), ier = 0 );
    }
  }
  else MPI_CHECK( MPI_Send(&mesh->point[1], mesh->np, MPI_POINT, parmesh->info.root, (parmesh->myrank+1)*MPI_ANALYS_TAG + parmesh->nprocs + 1, parmesh->comm), ier = 0 );

  // Metric
  nm_glob = np_glob * met->size;
  int displsm[parmesh->nprocs];
  for ( n = 0; n < parmesh->nprocs; n++ ) {
    list_np[n] *= met->size;
    displsm[n] = displsp[n] * met->size;
  }

  double* list_met = (double*) malloc(nm_glob * sizeof(double));
  if (parmesh->myrank == parmesh->info.root) {
    for ( n = 0; n < met->size*mesh->np; n++) list_met[n] = met->m[n+met->size];
    for ( n = 1; n < parmesh->nprocs; n++) {
      MPI_CHECK( MPI_Recv( &list_met[displsm[n]], list_np[n], MPI_DOUBLE, n, (n+1)*MPI_ANALYS_TAG + 2*(parmesh->nprocs + 1), parmesh->comm, &status), ier = 0 );
    }
  }
  else MPI_CHECK( MPI_Send(&met->m[met->size], met->size*mesh->np, MPI_DOUBLE, parmesh->info.root, (parmesh->myrank+1)*MPI_ANALYS_TAG + 2*(parmesh->nprocs + 1), parmesh->comm), ier = 0 );

  // Delete the MPI types
  PMMG2D_Free_MPI_meshDatatype(&MPI_POINT, &MPI_TRIANGLE, &MPI_TRIANGLE, 0);

  // 2) Fill the mesh on the root process
  if (parmesh->myrank != parmesh->info.root) return 1;

  PMMG2D_REALLOC(mesh, mesh->point, np_glob + 1, mesh->npmax + 1,
                 MMG5_Point,"Realloc point when merging the mesh", return 0);

  // Find the duplicate points in the mesh and delete them
  int* list_ref = (int*) malloc(np_glob * sizeof(int));
  int* index_points = (int*) malloc(np_glob * sizeof(int));
  for ( k = 0; k < np_glob; k++ ) list_ref[k] = 0;

  for ( k = 0; k < np_glob; k++ ) {
    if (!list_ref[list_points[k].ref]) {
      index_points[k] = k-offset+1;
      list_ref[list_points[k].ref] = index_points[k];
      mesh->point[index_points[k]] = list_points[k];
    }
    else {
      offset += 1;
      index_points[k] = list_ref[list_points[k].ref];
    }

  }

  mesh->np = np_glob - offset;
  met->np = np_glob - offset;

  PMMG2D_REALLOC(mesh, mesh->point, mesh->np + 1, np_glob + 1,
                 MMG5_Point, "Realloc point when merging the mesh (2)", return 0);

  // Build metric
  PMMG2D_REALLOC(mesh, met->m, met->size*(mesh->np+1), met->size * (mesh->npmax+1),
                 double,"Realloc metric when merging the mesh", return 0);

  for ( k = 0; k < np_glob; k++ ) {
    if (!list_ref[list_points[k].ref]) {
      for ( i = 0; i < 3; i++ ) met->m[met->size*index_points[k]+i] = list_met[met->size*k+i];
    }
  }

  free(list_points);
  free(list_met);
  free(list_ref);

  // Build triangles
  PMMG2D_REALLOC(mesh, mesh->tria, nt_glob + 1, mesh->ntmax + 1,
                 MMG5_Tria,"Realloc triangle when merging the mesh", return 0);

  mesh->nt = nt_glob;
  mesh->ntmax = mesh->nti = mesh->nt;
  mesh->xtmax = mesh->ntmax;

  // Renumber the vertex numbers in the triangle to be consistent with point numbers
  n = 1;
  offset = 0;
  for ( k = 0; k < nt_glob; k++ ) {

    if ( !MG_EOK(&list_triangles[k]) ) continue;

    mesh->tria[k+1] = list_triangles[k];

    if ( n < parmesh->nprocs ) {
      if ( k == displs[n] ) {
        offset = displsp[n];
        n++;
        while (n < parmesh->nprocs && list_np[n-1] == 0) {n++;}
      }
    }

    for ( i = 0; i < 3; i++ ) list_triangles[k].v[i] += offset;

  }

  for ( k = 1; k < mesh->nt+1; k++ ) {
    for ( i = 0; i < 3; i++ ) mesh->tria[k].v[i] = index_points[list_triangles[k-1].v[i]-1];
  }

  free(list_triangles);
  free(index_points);

  mesh->npi = mesh->np;
  mesh->npnil = 0;
  met->npi = mesh->np;
  mesh->npmax = mesh->np;
  met->npmax = met->np;

  // 3) Rebuild the boundary edges
  if (mesh->edge) {
    PMMG2D_DEL_MEM(mesh, mesh->edge, MMG5_Edge, "Delete edges when merging the mesh");
    PMMG2D_REALLOC(mesh, mesh->edge, 1, 0, MMG5_Edge, "Realloc edge when merging the mesh (1)", return 0);
    mesh->na = 0;
  }

  for ( k = 1; k <= mesh->nt; k++ ) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) continue;

    for ( i = 0; i < 3; i++ ) {

      if (pt->tag[i] & MG_BDY) {
        mesh->na++;
        PMMG2D_REALLOC(mesh, mesh->edge, mesh->na+1, mesh->na,
                       MMG5_Edge, "Realloc edge when merging the mesh (2)", return 0);
        i1 = MMG5_inxt2[i];
        i2 = MMG5_inxt2[i1];
        mesh->edge[mesh->na].a = pt->v[i1];
        mesh->edge[mesh->na].b = pt->v[i2];
        mesh->edge[mesh->na].tag = pt->tag[i];
        mesh->edge[mesh->na].ref = mesh->na;
        pt->edg[i] = mesh->na;
      }
    }
  }

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if fail, 1 if successc(on all procs)
 *
 *  merge all meshes and metrics to single mesh and metric in P0's parmesh
 */
int PMMG2D_merge_parmesh_initial( PMMG2D_pParMesh parmesh )
{
  int n, k, i, i1, i2, ier;
  MMG5_pMesh mesh, mesh_glob;
  MMG5_pSol met, met_glob;
  MMG5_pTria pt;
  MPI_Datatype   MPI_POINT, MPI_TRIANGLE;
  MPI_Status     status;
  int list_nt[parmesh->nprocs], list_np[parmesh->nprocs];
  int nt_glob = 0, np_glob = 0, nm_glob = 0;
  int offset = 0;

  mesh = parmesh->mesh;
  met = parmesh->met;

  // 1) Gather the mesh and metric information in root process
  for (n = 0; n < parmesh->nprocs; n++) {
    list_nt[n] = 0;
    list_np[n] = 0;
  }

  // Create specific types to MPI exchanges
  PMMG2D_create_MPI_Tria( &MPI_TRIANGLE );
  PMMG2D_create_MPI_Point( &MPI_POINT );

  // Triangles
  MPI_CHECK( MPI_Gather(&mesh->nt, 1, MPI_INT, list_nt, 1, MPI_INT, parmesh->info.root, parmesh->comm), ier = 0 );

  for ( n = 0; n < parmesh->nprocs; n++ ) nt_glob += list_nt[n];
  int displs[parmesh->nprocs];
  displs[0] = 0;
  for ( n = 1; n < parmesh->nprocs; n++ ) displs[n] = displs[n-1] + list_nt[n-1];

  // Problem when using MPI_Gather(v) with great meshes (probably due to MPI memory issues).
  // Send and Recv are therefore used so that MPI deals with smaller arrays each exchange.
  MMG5_pTria list_triangles = (MMG5_pTria) malloc(nt_glob * sizeof(MMG5_Tria));
  if (parmesh->myrank == parmesh->info.root) {
    for ( n = 0; n < mesh->nt; n++) list_triangles[n] = mesh->tria[n+1];
    for ( n = 1; n < parmesh->nprocs; n++) {
      MPI_CHECK( MPI_Recv( &list_triangles[displs[n]], list_nt[n], MPI_TRIANGLE, n, (n+1)*MPI_ANALYS_TAG, parmesh->comm, &status), ier = 0 );
    }
  }
  else MPI_CHECK( MPI_Send(&mesh->tria[1], mesh->nt, MPI_TRIANGLE, parmesh->info.root, (parmesh->myrank+1)*MPI_ANALYS_TAG, parmesh->comm), ier = 0 );

  // Points
  MPI_CHECK( MPI_Gather(&mesh->np, 1, MPI_INT, list_np, 1, MPI_INT, parmesh->info.root, parmesh->comm), ier = 0 );

  for ( n = 0; n < parmesh->nprocs; n++ ) np_glob += list_np[n];
  int displsp[parmesh->nprocs];
  displsp[0] = 0;
  for ( n = 1; n < parmesh->nprocs; n++ ) displsp[n] = displsp[n-1] + list_np[n-1];

  MMG5_pPoint list_points = (MMG5_pPoint) malloc(np_glob * sizeof(MMG5_Point));
  if (parmesh->myrank == parmesh->info.root) {
    for ( n = 0; n < mesh->np; n++) list_points[n] = mesh->point[n+1];
    for ( n = 1; n < parmesh->nprocs; n++) {
      MPI_CHECK( MPI_Recv( &list_points[displsp[n]], list_np[n], MPI_POINT, n, (n+1)*MPI_ANALYS_TAG + parmesh->nprocs + 1, parmesh->comm, &status), ier = 0 );
    }
  }
  else MPI_CHECK( MPI_Send(&mesh->point[1], mesh->np, MPI_POINT, parmesh->info.root, (parmesh->myrank+1)*MPI_ANALYS_TAG + parmesh->nprocs + 1, parmesh->comm), ier = 0 );

  // Metric
  nm_glob = np_glob * met->size;
  int displsm[parmesh->nprocs];
  for ( n = 0; n < parmesh->nprocs; n++ ) {
    list_np[n] *= met->size;
    displsm[n] = displsp[n] * met->size;
  }

  double* list_met = (double*) malloc(nm_glob * sizeof(double));
  if (parmesh->myrank == parmesh->info.root) {
    for ( n = 0; n < met->size*mesh->np; n++) list_met[n] = met->m[n+met->size];
    for ( n = 1; n < parmesh->nprocs; n++) {
      MPI_CHECK( MPI_Recv( &list_met[displsm[n]], list_np[n], MPI_DOUBLE, n, (n+1)*MPI_ANALYS_TAG + 2*(parmesh->nprocs + 1), parmesh->comm, &status), ier = 0 );
    }
  }
  else MPI_CHECK( MPI_Send(&met->m[met->size], met->size*mesh->np, MPI_DOUBLE, parmesh->info.root, (parmesh->myrank+1)*MPI_ANALYS_TAG + 2*(parmesh->nprocs + 1), parmesh->comm), ier = 0 );

  // Delete the MPI types
  PMMG2D_Free_MPI_meshDatatype(&MPI_POINT, &MPI_TRIANGLE, &MPI_TRIANGLE, 0);

  // 2) Fill the mesh on the root process
  if (parmesh->myrank != parmesh->info.root) return 1;

  if (parmesh->initial_mesh)  MMG2D_Free_all( MMG5_ARG_start,
                                              MMG5_ARG_ppMesh, &parmesh->initial_mesh,
                                              MMG5_ARG_ppMet,  &parmesh->initial_met,
                                              MMG5_ARG_end );

  parmesh->initial_mesh = NULL;
  parmesh->initial_met  = NULL;

  MMG2D_Init_mesh( MMG5_ARG_start,
                   MMG5_ARG_ppMesh, &parmesh->initial_mesh,
                   MMG5_ARG_ppMet, &parmesh->initial_met,
                   MMG5_ARG_end );

  // Set maximum memory
  parmesh->initial_mesh->memMax = parmesh->memGloMax;
  mesh_glob = parmesh->initial_mesh;
  mesh_glob->np  = 0;
  mesh_glob->nt  = 0;
  mesh_glob->npmax  = 0;
  mesh_glob->ntmax  = 0;
  met_glob = parmesh->initial_met;

  if ( !MMG2D_Set_solSize(parmesh->initial_mesh,parmesh->initial_met,MMG5_Vertex,0,parmesh->met->type) ) return 0;

  PMMG2D_REALLOC(mesh_glob, mesh_glob->point, np_glob + 1, mesh_glob->npmax + 1,
                 MMG5_Point,"Realloc point when merging the mesh", return 0);

  // Find the duplicate points in the mesh and delete them
  int* list_ref = (int*) malloc(np_glob * sizeof(int));
  int* index_points = (int*) malloc(np_glob * sizeof(int));
  for ( k = 0; k < np_glob; k++ ) list_ref[k] = 0;

  for ( k = 0; k < np_glob; k++ ) {

    if (!list_ref[list_points[k].ref]) {
      index_points[k] = k-offset+1;
      list_ref[list_points[k].ref] = index_points[k];
      mesh_glob->point[index_points[k]] = list_points[k];
    }
    else {
      offset += 1;
      index_points[k] = list_ref[list_points[k].ref];
    }

  }

  mesh_glob->np = np_glob - offset;
  met_glob->np = np_glob - offset;

  PMMG2D_REALLOC(mesh_glob, mesh_glob->point, mesh_glob->np + 1, np_glob + 1,
                 MMG5_Point, "Realloc point when merging the mesh (2)", return 0);

  // Build metric
  PMMG2D_REALLOC(mesh_glob, met_glob->m, met_glob->size*(mesh_glob->np+1), met_glob->size * (mesh_glob->npmax+1),
                 double,"Realloc metric when merging the mesh", return 0);

  for ( k = 0; k < np_glob; k++ ) {
    if (!list_ref[list_points[k].ref]) {
      for ( i = 0; i < 3; i++ ) met_glob->m[met_glob->size*index_points[k]+i] = list_met[met_glob->size*k+i];
    }
  }

  free(list_points);
  free(list_met);
  free(list_ref);

  // Build triangles
  PMMG2D_REALLOC(mesh_glob, mesh_glob->tria, nt_glob + 1, mesh_glob->ntmax + 1,
                 MMG5_Tria,"Realloc triangle when merging the mesh", return 0);

  mesh_glob->nt = nt_glob;
  mesh_glob->ntmax = mesh_glob->nti = mesh_glob->nt;
  mesh_glob->xtmax = mesh_glob->ntmax;

  // Renumber the vertex numbers in the triangle to be consistent with point numbers
  n = 1;
  offset = 0;
  for ( k = 0; k < nt_glob; k++ ) {

    if ( !MG_EOK(&list_triangles[k]) ) continue;

    mesh_glob->tria[k+1] = list_triangles[k];

    if ( n < parmesh->nprocs ) {
      if ( k == displs[n] ) {
        offset = displsp[n];
        n++;
        while (n < parmesh->nprocs && list_np[n-1] == 0) {n++;}
      }
    }

    for ( i = 0; i < 3; i++ ) list_triangles[k].v[i] += offset;

  }

  for ( k = 1; k < mesh_glob->nt+1; k++ ) {
    for ( i = 0; i < 3; i++ ) mesh_glob->tria[k].v[i] = index_points[list_triangles[k-1].v[i]-1];
  }

  mesh_glob->npi = mesh_glob->np;
  mesh_glob->npnil = 0;
  met_glob->npi = mesh_glob->np;
  mesh_glob->npmax = mesh_glob->np;
  met_glob->npmax = met_glob->np;

  free(list_triangles);
  free(index_points);

  return 1;
}
