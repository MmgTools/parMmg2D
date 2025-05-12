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
 * \file renumbermesh_pmmg2d.c
 * \brief Renumber the global node references over the processors.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */
#include "parmmg2d.h"

static int min(int array_[], int size)
{
  int min = array_[0];
  for (int i = 1; i < size; i++) min = fmin(min, array_[i]);
  return min;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 *
 * \return 1 if success, 0 otherwise.
 *
 * Renumber the vertices in the local meshes taking into account
 * the interface points which must have the same ref number
 * in the different processors.
 * Function for centralized mesh as input.
 */
int PMMG2D_renumber_mesh( PMMG2D_pParMesh parmesh )
{
  int i, n, j, k, size, ier = 1, ieresult;
  int* inverse_interface_node_list = (int*) malloc((parmesh->mesh->np+1) * sizeof(int));
  int rcounts[parmesh->nprocs], size_exchange[parmesh->nprocs];
  int** list_interface_nodes = (int**) malloc(parmesh->nprocs * sizeof(int*)); 
  for (i = 0; i < parmesh->nprocs; i++) list_interface_nodes[i] = (int*) malloc(sizeof(int));
  int** list_interface_local_index = (int**) malloc(parmesh->nprocs * sizeof(int*));
  for (i = 0; i < parmesh->nprocs; i++) list_interface_local_index[i] = (int*) malloc(sizeof(int));
  MPI_Status     status;

  MMG5_pMesh mesh;

  mesh = parmesh->mesh;

  // 1) Update the list of interface nodes (for the second step)
  PMMG2D_free_interface_nodes_list( parmesh );
  PMMG2D_fill_interface_nodes_list( parmesh );

  for (i = 1; i <= mesh->np; i++) inverse_interface_node_list[i] = 0;

  for (i = 0; i < parmesh->nip; i++) 
    inverse_interface_node_list[parmesh->list_interface_nodes[i].index[0]] = i;

  // 2) Count the number of nodes on each process
  // The interface nodes are added only on the process
  // with the smallest rank number
  n = 0;
  for (i = 1; i <= mesh->np; i++) {
    if (mesh->point[i].tag & MG_PARBDY) {
      size = parmesh->list_interface_nodes[inverse_interface_node_list[i]].nb;
      if (parmesh->myrank == min(parmesh->list_interface_nodes[inverse_interface_node_list[i]].proc, size)) n++;
    }
    else n++;
  }

  // 3) Gather this number of nodes on every process
  MPI_Allgather(&n, 1, MPI_INT, rcounts, 1, MPI_INT, parmesh->comm );

  // 4) Renumber the vertices
  for (i = 0; i < parmesh->nprocs; i++) size_exchange[i] = 0;
  n = 0;
  for (i = 0; i < parmesh->myrank; i++) n += rcounts[i];
  
  j = 1;
  for (i = 1; i <= mesh->np; i++) {
    if (mesh->point[i].tag & MG_PARBDY) {
      size = parmesh->list_interface_nodes[inverse_interface_node_list[i]].nb;
      if (parmesh->myrank == min(parmesh->list_interface_nodes[inverse_interface_node_list[i]].proc, size)) {
        mesh->point[i].ref = n + j;
        for (k = 1; k < size; k++) {
          int p = parmesh->list_interface_nodes[inverse_interface_node_list[i]].proc[k];
          list_interface_nodes[p][size_exchange[p]] = n+j;
          list_interface_nodes[p] = (int*) realloc(list_interface_nodes[p], (size_exchange[p]+2) * sizeof(int));
          list_interface_local_index[p][size_exchange[p]] = 
            parmesh->list_interface_nodes[inverse_interface_node_list[i]].index[k];
          list_interface_local_index[p] = (int*) realloc(list_interface_local_index[p], (size_exchange[p]+2) * sizeof(int));
          size_exchange[p]++;
        }
        j++;
      }
    }
    else {
      mesh->point[i].ref = n + j;
      j++;
    }
  }

  free(inverse_interface_node_list);

  // 5) Exchange the ref number of PARBDY nodes to neighbours
  for (i = 0; i < parmesh->nprocs; i++) {

    if (i == parmesh->myrank) continue;

    MPI_CHECK( MPI_Send(&size_exchange[i], 1, MPI_INT, i, MPI_ANALYS_TAG + i,
                         parmesh->comm), ier = 0 );

    MPI_CHECK( MPI_Recv(&size, 1, MPI_INT, i, MPI_ANALYS_TAG + parmesh->myrank,
                         parmesh->comm, &status), ier = 0 );

    if (size_exchange[i] > 0) {
      MPI_CHECK( MPI_Send(&list_interface_nodes[i][0], size_exchange[i], MPI_INT, i, MPI_ANALYS_TAG + i + parmesh->nprocs,
                           parmesh->comm), ier = 0 );
      MPI_CHECK( MPI_Send(&list_interface_local_index[i][0], size_exchange[i], MPI_INT, i, MPI_ANALYS_TAG + i + 2*parmesh->nprocs,
                           parmesh->comm), ier = 0 );
    }

    if (size > 0) {
      int* list_interface_nodess = (int*) malloc(size * sizeof(int));
      MPI_CHECK( MPI_Recv(&list_interface_nodess[0], size, MPI_INT, i, MPI_ANALYS_TAG + parmesh->myrank + parmesh->nprocs,
                           parmesh->comm, &status), ier = 0 );
      int* list_interface_local_indices = (int*) malloc(size * sizeof(int));
      MPI_CHECK( MPI_Recv(&list_interface_local_indices[0], size, MPI_INT, i, MPI_ANALYS_TAG + parmesh->myrank + 2*parmesh->nprocs,
                           parmesh->comm, &status), ier = 0 );

      for (j = 0; j < size; j++) mesh->point[list_interface_local_indices[j]].ref = list_interface_nodess[j];

      free(list_interface_nodess);
      free(list_interface_local_indices);
    }
  }

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

  for (i = 0; i < parmesh->nprocs; i++) {
    free(list_interface_nodes[i]);
    free(list_interface_local_index[i]);
  }

  free(list_interface_nodes);
  free(list_interface_local_index);

  return ieresult;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 *
 * \return 1 if success, 0 otherwise.
 *
 * Renumber the vertices in the local meshes taking into account
 * the interface points which must have the same ref number
 * in the different processors.
 * Function for distributed meshes as input.
 */
int PMMG2D_renumber_distributed_mesh( PMMG2D_pParMesh parmesh, PMMG2D_pint_nodes list_interface_nodes, int size_interface )
{
  int i, n, j, k, size, ier, ieresult;
  int* inverse_interface_node_list = (int*) malloc((parmesh->mesh->np+1) * sizeof(int));
  int rcounts[parmesh->nprocs], size_exchange[parmesh->nprocs];
  int** list_interface_local_nodes = (int**) malloc(parmesh->nprocs * sizeof(int*));
  for (i = 0; i < parmesh->nprocs; i++) list_interface_local_nodes[i] = (int*) malloc(sizeof(int));
  int** list_interface_local_index = (int**) malloc(parmesh->nprocs * sizeof(int*));
  for (i = 0; i < parmesh->nprocs; i++) list_interface_local_index[i] = (int*) malloc(sizeof(int));

  MPI_Status     status;

  MMG5_pMesh mesh;

  mesh = parmesh->mesh;

  for (i = 1; i <= mesh->np; i++) inverse_interface_node_list[i] = 0;

  for (i = 0; i < size_interface; i++)
    inverse_interface_node_list[list_interface_nodes[i].index[0]] = i;

  // 2) Count the number of nodes on each process
  // The interface nodes are added only on the process
  // with the smallest rank number
  n = 0;
  for (i = 1; i <= mesh->np; i++) {
    if (mesh->point[i].tag & MG_PARBDY) {
      size = list_interface_nodes[inverse_interface_node_list[i]].nb;
      if (parmesh->myrank == min(list_interface_nodes[inverse_interface_node_list[i]].proc, size)) n++;
    }
    else n++;
  }

  // 3) Gather this number of nodes on every process
  MPI_Allgather(&n, 1, MPI_INT, rcounts, 1, MPI_INT, parmesh->comm );

  // 4) Renumber the vertices
  for (i = 0; i < parmesh->nprocs; i++) size_exchange[i] = 0;
  n = 0;
  for (i = 0; i < parmesh->myrank; i++) n += rcounts[i];

  j = 1;
  for (i = 1; i <= mesh->np; i++) {
    if (mesh->point[i].tag & MG_PARBDY) {
      size = list_interface_nodes[inverse_interface_node_list[i]].nb;
      if (parmesh->myrank == min(list_interface_nodes[inverse_interface_node_list[i]].proc, size)) {
        mesh->point[i].ref = n + j;
        for (k = 1; k < size; k++) {
          int p = list_interface_nodes[inverse_interface_node_list[i]].proc[k];
          list_interface_local_nodes[p][size_exchange[p]] = n+j;
          list_interface_local_nodes[p] = (int*) realloc(list_interface_local_nodes[p], (size_exchange[p]+2) * sizeof(int));
          list_interface_local_index[p][size_exchange[p]] =
            list_interface_nodes[inverse_interface_node_list[i]].index[k];
          list_interface_local_index[p] = (int*) realloc(list_interface_local_index[p], (size_exchange[p]+2) * sizeof(int));
          size_exchange[p]++;
        }
        j++;
      }
    }
    else {
      mesh->point[i].ref = n + j;
      j++;
    }
  }

  free(inverse_interface_node_list);

  // 5) Exchange the ref number of PARBDY nodes to neighbours
  for (i = 0; i < parmesh->nprocs; i++) {

    if (i == parmesh->myrank) continue;

    MPI_CHECK( MPI_Send(&size_exchange[i], 1, MPI_INT, i, MPI_ANALYS_TAG + i,
                         parmesh->comm), ier = 0 );

    MPI_CHECK( MPI_Recv(&size, 1, MPI_INT, i, MPI_ANALYS_TAG + parmesh->myrank,
                         parmesh->comm, &status), ier = 0 );

    if (size_exchange[i] > 0) {
      MPI_CHECK( MPI_Send(&list_interface_local_nodes[i][0], size_exchange[i], MPI_INT, i, MPI_ANALYS_TAG + i + parmesh->nprocs,
                           parmesh->comm), ier = 0 );
      MPI_CHECK( MPI_Send(&list_interface_local_index[i][0], size_exchange[i], MPI_INT, i, MPI_ANALYS_TAG + i + 2*parmesh->nprocs,
                           parmesh->comm), ier = 0 );
    }

    if (size > 0) {
      int* list_interface_nodes_ref = (int*) malloc(size * sizeof(int));
      MPI_CHECK( MPI_Recv(&list_interface_nodes_ref[0], size, MPI_INT, i, MPI_ANALYS_TAG + parmesh->myrank + parmesh->nprocs,
                           parmesh->comm, &status), ier = 0 );
      int* list_interface_local_index = (int*) malloc(size * sizeof(int));
      MPI_CHECK( MPI_Recv(&list_interface_local_index[0], size, MPI_INT, i, MPI_ANALYS_TAG + parmesh->myrank + 2*parmesh->nprocs,
                           parmesh->comm, &status), ier = 0 );

      for (j = 0; j < size; j++) mesh->point[list_interface_local_index[j]].ref = list_interface_nodes_ref[j];

      free(list_interface_nodes_ref);
      free(list_interface_local_index);
    }
  }

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

  for (i = 0; i < parmesh->nprocs; i++) {
    free(list_interface_local_nodes[i]);
    free(list_interface_local_index[i]);
  }

  free(list_interface_local_nodes);
  free(list_interface_local_index);

  return ieresult;
}

