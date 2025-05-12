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
 * \file PT_Scotch_pmmg2d.c
 * \brief Partition mesh using Scotch and PT-Scotch
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "parmmg2d.h"

#ifdef USE_SCOTCH
#include <scotch.h>
#endif

#ifdef USE_PTSCOTCH
#include <ptscotch.h>
#endif

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param nb_neighbours list of the numbers of neighbour triangles for each triangle.
 * \param list_neighbours list of neighbour triangles (global reference number).
 * \return 0 if it fails, 1 otherwise.
 *
 * Find the neighbour triangles of each triangle (including "ghost" triangles).
 *
 */
int PMMG2D_find_neighbour_triangles( PMMG2D_pParMesh parmesh, int* nb_neighbours, int **list_neighbours )
{
  MMG5_pMesh mesh;
  MMG5_pTria pt, pt_neigh;
  int ier = 1, ieresult;
  MPI_Status status;

  mesh = parmesh->mesh;

  // Change the reference number of the triangles to be consecutive and globally unique
  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if (!MG_EOK(pt)) continue;

    pt->ref = k-1;
  }

  int rcounts[parmesh->nprocs];
  MPI_Allgather(&mesh->nt, 1, MPI_INT, rcounts, 1, MPI_INT, parmesh->comm);

  int n = 0;
  for (int k = 0; k < parmesh->myrank; k++) n += rcounts[k];

  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];
    if (!MG_EOK(pt)) continue;
    pt->ref += n;
  }

  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) )  continue;

    for(int i = 1; i < 4; i++ ) { // Loop over the adjacent triangles to triangle k

      int k1 = mesh->adja[3*(k-1)+i];

      if( !k1 ) continue;

      k1 /= 3;
      pt_neigh = &mesh->tria[k1];

      list_neighbours[k-1][nb_neighbours[k-1]++] = pt_neigh->ref;

    }
  }

  // Store the interface nodes
  // Greatest reference number
  int max_ref = 0;
  for (int i = 1; i <= mesh->np; i++) {
    if (mesh->point[i].ref > max_ref) max_ref = mesh->point[i].ref;
  }
  int max_ref_global;
  MPI_Allreduce( &max_ref, &max_ref_global, 1, MPI_INT, MPI_MAX, parmesh->comm );

  int* interface_nodes = (int*) malloc(max_ref_global * sizeof(int));
  for (int k = 0; k < parmesh->nip; k++) {
    interface_nodes[mesh->point[(&parmesh->list_interface_nodes[k])->index[0]].ref] = k;
  }

  // Build the list of the neighbouring triangles at the parallel interfaces
  int size[parmesh->nprocs];
  for (int k = 0; k < parmesh->nprocs; k++) size[k] = 0;

  int **list_neighbours_interface = malloc(parmesh->nprocs*sizeof(int*));
  for (int k = 0; k < parmesh->nprocs; k++) list_neighbours_interface[k] = malloc(1*sizeof(int));

  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) )  continue;

    for (int i = 0; i < 3; i++) {

      // If the edge is not on an parallel interface
      if (!(pt->tag[i] & MG_PARBDY)) continue;

      // Processors containing the first point
      int k1 = interface_nodes[mesh->point[pt->v[MMG5_inxt2[i]]].ref];
      int list_proc1[(&parmesh->list_interface_nodes[k1])->nb-1];
      for (int j = 1; j < (&parmesh->list_interface_nodes[k1])->nb; j++) {
        list_proc1[j-1] = (&parmesh->list_interface_nodes[k1])->proc[j];
      }

      // Processors containing the second point
      int k2 = interface_nodes[mesh->point[pt->v[MMG5_iprv2[i]]].ref];
      int list_proc2[(&parmesh->list_interface_nodes[k2])->nb-1];
      for (int j = 1; j < (&parmesh->list_interface_nodes[k2])->nb; j++) {
        list_proc2[j-1] = (&parmesh->list_interface_nodes[k2])->proc[j];
      }

      // Find the common processor of both points
      int should_break = 0;
      int neighbour_proc;
      for (int j = 0; j < (&parmesh->list_interface_nodes[k1])->nb-1; j++) {
        for (int l = 0; l < (&parmesh->list_interface_nodes[k2])->nb-1; l++) {
          if (list_proc1[j] == list_proc2[l]) {
            neighbour_proc = list_proc1[j];
            should_break = 1;
            break;
          }
        }
        if (should_break) break;
      }

      // Add the triangle and the point references to the exchange list associate with processor "neighbour_proc"
      list_neighbours_interface[neighbour_proc] = (int*)realloc(list_neighbours_interface[neighbour_proc], (size[neighbour_proc] + 3)*sizeof(int));
      list_neighbours_interface[parmesh->myrank] = (int*)realloc(list_neighbours_interface[parmesh->myrank], (size[parmesh->myrank] + 3)*sizeof(int));
      list_neighbours_interface[neighbour_proc][size[neighbour_proc]++] = pt->ref;
      list_neighbours_interface[parmesh->myrank][size[parmesh->myrank]++] = k-1;

      if (mesh->point[pt->v[MMG5_iprv2[i]]].ref < mesh->point[pt->v[MMG5_inxt2[i]]].ref) { // Ascending order to simplify the search below
        list_neighbours_interface[neighbour_proc][size[neighbour_proc]++] = mesh->point[pt->v[MMG5_iprv2[i]]].ref;
        list_neighbours_interface[neighbour_proc][size[neighbour_proc]++] = mesh->point[pt->v[MMG5_inxt2[i]]].ref;
        list_neighbours_interface[parmesh->myrank][size[parmesh->myrank]++] = mesh->point[pt->v[MMG5_iprv2[i]]].ref;
        list_neighbours_interface[parmesh->myrank][size[parmesh->myrank]++] = mesh->point[pt->v[MMG5_inxt2[i]]].ref;
      }
      else {
        list_neighbours_interface[neighbour_proc][size[neighbour_proc]++] = mesh->point[pt->v[MMG5_inxt2[i]]].ref;
        list_neighbours_interface[neighbour_proc][size[neighbour_proc]++] = mesh->point[pt->v[MMG5_iprv2[i]]].ref;
        list_neighbours_interface[parmesh->myrank][size[parmesh->myrank]++] = mesh->point[pt->v[MMG5_inxt2[i]]].ref;
        list_neighbours_interface[parmesh->myrank][size[parmesh->myrank]++] = mesh->point[pt->v[MMG5_iprv2[i]]].ref;
      }

    }
  }

  free(interface_nodes);

  // Exchange the list of neighbouring triangles
  int size_list;
  for (int i = 0; i < parmesh->nprocs; i++) {

    if (i == parmesh->myrank) continue;

    MPI_CHECK( MPI_Send(&size[i], 1, MPI_INT, i, MPI_ANALYS_TAG + i, parmesh->comm), ier = 0 );

    MPI_CHECK( MPI_Recv(&size_list, 1, MPI_INT, i, MPI_ANALYS_TAG + parmesh->myrank, parmesh->comm, &status), ier = 0 );

    if (size[i] > 0) MPI_CHECK( MPI_Send(&list_neighbours_interface[i][0], size[i], MPI_INT, i,
                                         MPI_ANALYS_TAG + i + parmesh->nprocs, parmesh->comm), ier = 0 );

    if (size_list > 0) {
      int* list_interface_triangles = (int*) malloc(size_list * sizeof(int));
      MPI_CHECK( MPI_Recv(&list_interface_triangles[0], size_list, MPI_INT, i, MPI_ANALYS_TAG + parmesh->myrank + parmesh->nprocs,
                           parmesh->comm, &status), ier = 0 );

      // Loop over all the triangles received from proc i and add them
      int k = 0;
      while (k < size_list) {
        int n_triangle = list_interface_triangles[k++];
        int p1 = list_interface_triangles[k++];
        int p2 = list_interface_triangles[k++];

        for (int j = 0; j < size[parmesh->myrank]/3; j++) {

          // Find the neighbouring triangle which also contains the two interface nodes
          if (list_neighbours_interface[parmesh->myrank][3*j+1] != p1) continue;
          if (list_neighbours_interface[parmesh->myrank][3*j+2] != p2) continue;

          int ref = list_neighbours_interface[parmesh->myrank][3*j]; // local index
          list_neighbours[ref] = (int*)realloc(list_neighbours[ref], (nb_neighbours[ref]+1)*sizeof(int));
          list_neighbours[ref][nb_neighbours[ref]++] = n_triangle; // global index of the neighbour
          break;
        }

      }
      free(list_interface_triangles);
    }
  }

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

  for (int k = 0; k < parmesh->nprocs; k++) free(list_neighbours_interface[k]);
  free(list_neighbours_interface);

  return ieresult;
}

#ifdef USE_PTSCOTCH
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param partTab list of partitions of each triangle
 * \return 0 if it fails, 1 otherwise. 
 *
 * Modify the local mesh to fit with the new domain decomposition
 *
 */
int PMMG2D_PTScotch_exchange( PMMG2D_pParMesh parmesh, SCOTCH_Num *partTab )
{
  int list_index_pt[parmesh->nprocs];
  int list_index_pp[parmesh->nprocs];
  
  MMG5_pTria pt_global = NULL;
  MMG5_pPoint pp_global = NULL;
  double* mm_global = NULL;
  MMG5_pTria pt;
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;

  MMG5_pTria* list_pt = (MMG5_pTria*) malloc(parmesh->nprocs * sizeof(MMG5_pTria));
  for (int i = 0; i < parmesh->nprocs; i++) list_pt[i] = (MMG5_pTria) malloc(mesh->nt * sizeof(MMG5_Tria));
  MMG5_pPoint* list_pp = (MMG5_pPoint*) malloc(parmesh->nprocs * sizeof(MMG5_pPoint));
  for (int i = 0; i < parmesh->nprocs; i++) list_pp[i] = (MMG5_pPoint) malloc(3*mesh->nt * sizeof(MMG5_Point));
  double** list_mm = (double**) malloc(parmesh->nprocs * sizeof(double*));
  for (int i = 0; i < parmesh->nprocs; i++) list_mm[i] = (double*) malloc(parmesh->met->size*3*mesh->nt * sizeof(double));
  int* list_deleted_triangles = (int*) malloc((mesh->nt+1) * sizeof(int));

  int ier, ieresult;
  int npt_global, npp_global;

  for (int i = 0; i < parmesh->nprocs; i++) {
    list_index_pt[i] = 0;
    list_index_pp[i] = 0;
  }
  for (int i = 0; i < mesh->nt; i++) list_deleted_triangles[i] = 0;

  // 1) Find the triangles to exchange
  int n = 0;
  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    if (partTab[k-1] == parmesh->myrank) continue;

    int send_proc = partTab[k-1];

    list_deleted_triangles[n++] = k;

    list_pt[send_proc][list_index_pt[send_proc]++] = mesh->tria[k]; 

    for (int i = 0; i <= 2; i++) {
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

  // 2) Exchange the triangles, points and metrics
  ier = PMMG2D_exchange( parmesh, mesh->nt, list_index_pt, list_pt, list_index_pp, list_pp, list_mm,
                         &pt_global, &pp_global, &mm_global, &npt_global, &npp_global );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error in the mpi exchanges during the PT-Scotch decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  for (int i = 0; i < parmesh->nprocs; i++) {
    free(list_pt[i]);
    free(list_pp[i]);
    free(list_mm[i]);
  }
  free(list_pt);
  free(list_pp);
  free(list_mm);

  // 3) Remove the sent triangles from the local mesh
  int j = 0;
  for (int k = 1; k <= mesh->nt; k++ ) {
    if (list_deleted_triangles[j] == k) {
      j += 1;
      continue;
    }
    mesh->tria[k-j] = mesh->tria[k];
  }

  free(list_deleted_triangles);

  mesh->nt -= j;
  int nt_tmp = mesh->nt - j; // temporary number of triangles, used to check the contiguity (removing twice the triangles)

  PMMG2D_REALLOC(mesh, mesh->tria, mesh->nt+1, mesh->ntmax+1,
                 MMG5_Tria,"PT-Scotch decomposition realloc triangles",return 0);
  PMMG2D_REALLOC(mesh, mesh->point, mesh->np+1, mesh->npmax+1,
                 MMG5_Point,"PT-Scotch decomposition realloc points", return 0);
  PMMG2D_REALLOC(mesh, parmesh->met->m, parmesh->met->size*(mesh->np+1), parmesh->met->size*(mesh->npmax+1),
                 double, "PT-Scotch decomposition realloc metrics", return 0);

  // 4) Remove the points that do not belong to the mesh anymore
  ier = PMMG2D_remove_points( parmesh );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when removing the points during the PT-Scotch decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 5) Add the triangles to their corresponding new process
  ier = PMMG2D_add_triangles( parmesh, pt_global, pp_global, mm_global, npt_global, npp_global);
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when adding the new triangles during the PT-Scotch decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  PMMG2D_update( parmesh );

  // 6) Check the contiguity and correct the isolated triangles
  ier = PMMG2D_check_contiguity( parmesh, nt_tmp );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when checking the contiguity during the PT-Scotch decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  PMMG2D_update( parmesh );

  // 7) Untag the interface parallel points and edges
  PMMG2D_untag_parallel( parmesh );

  // 8) Mark the interface nodes
  ier = PMMG2D_mark_parallel_interface_nodes( parmesh );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when marking the interface nodes during the contiguity check in the PT-Scotch decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 9) Update the list of interface nodes
  PMMG2D_free_interface_nodes_list( parmesh );
  PMMG2D_fill_interface_nodes_list( parmesh );

  // Free the pointers
  PMMG2D_DEL_MEM(mesh, pt_global, MMG5_Tria, "Delete pt_global in the PT-Scotch decomposition");
  PMMG2D_DEL_MEM(mesh, pp_global, MMG5_Point, "Delete pp_global in the PT-Scotch decomposition");
  PMMG2D_DEL_MEM(mesh, mm_global, double, "Delete mm_global in the PT-Scotch decomposition");

  // 10) Rebuild the boundary edges
  if (mesh->edge) {
    PMMG2D_DEL_MEM(mesh, mesh->edge, MMG5_Edge, "Delete edges when decomposing the mesh using PT-Scotch");
    PMMG2D_REALLOC(mesh, mesh->edge, 1, 0, MMG5_Edge, "Realloc edge when decomposing the mesh using PT-Scotch (1)", return 0);
    mesh->na = 0;
  }

  for (int k = 1; k <= mesh->nt; k++ ) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) continue;

    for (int i = 0; i < 3; i++ ) {

      if (pt->tag[i] & MG_BDY) {
        mesh->na++;
        PMMG2D_REALLOC(mesh, mesh->edge, mesh->na+1, mesh->na,
                       MMG5_Edge, "Realloc edge when decomposing the mesh using PT-Scotch (2)", return 0);
        int i1 = MMG5_inxt2[i];
        int i2 = MMG5_inxt2[i1];
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
 * \return 0 if it fails, 1 otherwise. 
 * 
 * Decompose the mesh in parallel using PT-Scotch to balance the loads.
 *
 */
int PMMG2D_PTScotch_decomposition( PMMG2D_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  SCOTCH_Dgraph dgraph;
  SCOTCH_Num numVertices, numParts, *vertTab, *edgeTab, *partTab;
  SCOTCH_Num baseval = 0;
  int ier, ieresult;
  MPI_Status status;

  mesh = parmesh->mesh;

  numParts = parmesh->nprocs;

  // List of the triangles that are the neighbours
  int* nb_neighbours = malloc(mesh->nt*sizeof(int));
  int NB_TRIANGLES = 3; // Maximum number of neighbouring triangles associated to each triangle
  int **list_neighbours = malloc(mesh->nt*sizeof(int*));
  for (int i = 0; i < mesh->nt; i++) {
    nb_neighbours[i] = 0;
    list_neighbours[i] = malloc(NB_TRIANGLES*sizeof(int));
    for (int j = 0; j < NB_TRIANGLES; j++) list_neighbours[i][j] = 0;
  }

  if (!PMMG2D_find_neighbour_triangles(parmesh, nb_neighbours, list_neighbours )) {
    fprintf(stderr, "Error when finding the neighbouring triangles in PT-Scotch decomposition.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  numVertices = mesh->nt;
  vertTab = malloc((numVertices + 1) * sizeof(SCOTCH_Num));
  partTab = malloc(numVertices * sizeof(SCOTCH_Num));

  // Fill the points and edges structures
  int numEdges = 0;
  for (int k = 0; k < mesh->nt; k++) numEdges += nb_neighbours[k];
  edgeTab = malloc(numEdges * sizeof(SCOTCH_Num));

  vertTab[0] = 0;
  for (int k = 0; k < mesh->nt; k++) {
    vertTab[k+1] = vertTab[k] + nb_neighbours[k];
    for (int i = 0; i < nb_neighbours[k]; i++) {
      edgeTab[vertTab[k]+i] = list_neighbours[k][i];
    }
  }

  // Initialisation of the distributed graph
  if (SCOTCH_dgraphInit(&dgraph, parmesh->comm)) {
    fprintf(stderr, "Error when initialising the distributed graph.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Build the distributed graph
  if (SCOTCH_dgraphBuild(&dgraph, baseval,
                         numVertices, numVertices, vertTab, vertTab + 1,
                         NULL, NULL, numEdges, numEdges, edgeTab, NULL, NULL)) {
    fprintf(stderr, "Error when building the distributed graph.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Check the graph
  if (SCOTCH_dgraphCheck(&dgraph)) {
    fprintf(stderr, "Error in the distributed graph.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Definition of the strategy
  SCOTCH_Strat strat;
  SCOTCH_stratInit(&strat);
  if (SCOTCH_strategy == 1) {
    // For HPC, minimize the exchanges between process. 
    // The third argument is the number of attempts to optimize the partitioning. 
    // The last argument is the tolerance of imbalance in %.
    SCOTCH_stratDgraphMapBuild(&strat, SCOTCH_STRATQUALITY, parmesh->nprocs, parmesh->nprocs, PMMG2D_LOAD_IMBALANCE);
  }
  else if (SCOTCH_strategy == 2) {
    // To enforce a good balance between processors
    SCOTCH_stratDgraphMapBuild(&strat, SCOTCH_STRATBALANCE, parmesh->nprocs, parmesh->nprocs, PMMG2D_LOAD_IMBALANCE);
  }
  else if (SCOTCH_strategy == 3) { 
    SCOTCH_stratDgraphMapBuild(&strat, SCOTCH_STRATSPEED, parmesh->nprocs, parmesh->nprocs, PMMG2D_LOAD_IMBALANCE);
  }

  if (SCOTCH_dgraphPart(&dgraph, numParts, &strat, partTab)) {
    fprintf(stderr, "Error when partitioning the graph.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  SCOTCH_stratExit(&strat);
  SCOTCH_dgraphFree(&dgraph);

  for (int i = 0; i < mesh->nt; i++) free(list_neighbours[i]);
  free(list_neighbours);
  free(nb_neighbours);

  // Exchange the triangles to each new process
  ier = PMMG2D_PTScotch_exchange(parmesh,partTab);

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

  free(vertTab);
  free(edgeTab);
  free(partTab);

  return ieresult;
}
#endif

#ifdef USE_SCOTCH
/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * Decompose the mesh sequentially using Scotch to balance the loads.
 *
 */
int PMMG2D_Scotch_decomposition_root( PMMG2D_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  MMG5_pTria pt;
  SCOTCH_Graph graph;
  SCOTCH_Num numVertices, numParts, *vertTab, *edgeTab, *partTab;
  SCOTCH_Num baseval = 0;
  int ier, ieresult;
  MPI_Status status;

  mesh = parmesh->mesh;

  numParts = parmesh->nprocs;

  if ( (!mesh->adja) && (1 != MMG2D_hashTria( mesh )) ) {
    fprintf( stderr,"  ## PMMG2D Hashing problem in Scotch decomposition.\n" );
    return 0;
  }

  // List of the triangles that are the neighbours
  int* nb_neighbours = (int*) malloc(mesh->nt * sizeof(int));
  int NB_TRIANGLES = 3; // Maximum number of neighbouring triangles associated to each triangle
  int **list_neighbours = malloc(mesh->nt*sizeof(int*));
  for (int i = 0; i < mesh->nt; i++) {
    nb_neighbours[i] = 0;
    list_neighbours[i] = malloc(NB_TRIANGLES*sizeof(int));
    for (int j = 0; j < NB_TRIANGLES; j++) list_neighbours[i][j] = 0;
  }

  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) )  continue;

    for(int i = 1; i < 4; i++ ) { // Loop over the adjacent triangles to triangle k

      int k1 = mesh->adja[3*(k-1)+i];

      if( !k1 ) continue;

      k1 /= 3;

      list_neighbours[k-1][nb_neighbours[k-1]++] = k1-1;

    }
  }

  int numEdges = 0;
  for (int k = 0; k < mesh->nt; k++) numEdges += nb_neighbours[k];
  numVertices = mesh->nt;
  vertTab = malloc((numVertices + 1) * sizeof(SCOTCH_Num));
  partTab = malloc(numVertices * sizeof(SCOTCH_Num));
  edgeTab = malloc(numEdges * sizeof(SCOTCH_Num));

  // Fill the points and edges structures
  vertTab[0] = 0;
  for (int k = 0; k < mesh->nt; k++) {

    vertTab[k+1] = vertTab[k] + nb_neighbours[k];

    for (int i = 0; i < nb_neighbours[k]; i++) {
      edgeTab[vertTab[k]+i] = list_neighbours[k][i];
    }
  }

  free(nb_neighbours);

  // Initialisation of the distributed graph
  if (SCOTCH_graphInit(&graph)) {
    fprintf(stderr, "Error when initialising the centralized graph.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Build the distributed graph
  if (SCOTCH_graphBuild(&graph, baseval,
                         numVertices, vertTab, vertTab + 1,
                         NULL, NULL, numEdges, edgeTab, NULL)) {
    fprintf(stderr, "Error when building the centralized graph.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Check the graph
  if (SCOTCH_graphCheck(&graph)) {
    fprintf(stderr, "Error in the centralized graph.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Definition of the strategy
  SCOTCH_Strat strat;
  SCOTCH_stratInit(&strat);

  if (SCOTCH_strategy == 1) {
    // For HPC, minimize the exchanges between process. 
    // The third argument is the number of attempts to optimize the partitioning. 
    // The last argument is the tolerance of imbalance in %.
    SCOTCH_stratGraphMapBuild(&strat, SCOTCH_STRATQUALITY, parmesh->nprocs, PMMG2D_LOAD_IMBALANCE);
  }
  else if (SCOTCH_strategy == 2) {
    // To enforce a good balance between processors
    SCOTCH_stratGraphMapBuild(&strat, SCOTCH_STRATBALANCE, parmesh->nprocs, PMMG2D_LOAD_IMBALANCE);
  }
  else if (SCOTCH_strategy == 3) { 
    SCOTCH_stratGraphMapBuild(&strat, SCOTCH_STRATSPEED, parmesh->nprocs, PMMG2D_LOAD_IMBALANCE);
  }

  if (SCOTCH_graphPart(&graph, numParts, &strat, partTab)) {
    fprintf(stderr, "Error when partitioning the graph.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  SCOTCH_stratExit(&strat);
  SCOTCH_graphFree(&graph);

  for (int i = 0; i < mesh->nt; i++) free(list_neighbours[i]);
  free(list_neighbours);

  // Store the partition number of each triangle in base
  for (int k = 1; k <= mesh->nt; k++) mesh->tria[k].base = partTab[k-1];

  free(vertTab);
  free(edgeTab);
  free(partTab);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if it fails, 1 otherwise. 
 *
 * Decompose the mesh sequentially using Scotch and spread it over processors.
 *
 */
int PMMG2D_Scotch_decomposition( PMMG2D_pParMesh parmesh )
{
  int ier, ieresult;
  MMG5_pMesh mesh;
  MMG5_pTria pt;

  mesh = parmesh->mesh;

  // Decompose the mesh
  if ( !parmesh->myrank ) ier = PMMG2D_Scotch_decomposition_root(parmesh);

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error in the Scotch decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  if (mesh->edge) {
    PMMG2D_DEL_MEM(mesh, mesh->edge, MMG5_Edge, "Delete edges when decomposing the mesh using Scotch");
    PMMG2D_REALLOC(mesh, mesh->edge, 1, 0, MMG5_Edge, "Realloc edge when decomposing the mesh using Scotch (1)", return 0);
    mesh->na = 0;
  }

  // Send the local meshes from root process to the corresponding procs
  ier = PMMG2D_Metis_exchange(parmesh); // Same function as with Metis

  // Rebuild the boundary edges
  for (int k = 1; k <= mesh->nt; k++ ) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) continue;

    for (int i = 0; i < 3; i++ ) {

      if (pt->tag[i] & MG_BDY) {
        mesh->na++;
        PMMG2D_REALLOC(mesh, mesh->edge, mesh->na+1, mesh->na,
                       MMG5_Edge, "Realloc edge when decomposing the mesh using Scotch (2)", return 0);
        int i1 = MMG5_inxt2[i];
        int i2 = MMG5_inxt2[i1];
        mesh->edge[mesh->na].a = pt->v[i1];
        mesh->edge[mesh->na].b = pt->v[i2];
        mesh->edge[mesh->na].tag = pt->tag[i];
        mesh->edge[mesh->na].ref = mesh->na;
        pt->edg[i] = mesh->na;
      }
    }
  }
  
  return ier;
}
#endif
