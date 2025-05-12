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
 * \file parallel_interface.c
 * \brief Handle the parallel interface between the processors.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 * \todo doxygen documentation.
 */
#include "parmmg2d.h"

// Structure used in PMMG2D_fill_interface_nodes_list
typedef struct {
  int proc;
  int global_index;
  int local_index;
} InterfaceNode;

// Comparison functions used in PMMG2D_fill_interface_nodes_list
int compare_interface_nodes(const void *a, const void *b) {
  const InterfaceNode *elementA = (const InterfaceNode *)a;
  const InterfaceNode *elementB = (const InterfaceNode *)b;
  return (elementA->global_index > elementB->global_index) - (elementA->global_index < elementB->global_index);
}

int compare_int_nodes(const void *a, const void *b) {
  const PMMG2D_int_nodes *elementA = (const PMMG2D_int_nodes *)a;
  const PMMG2D_int_nodes *elementB = (const PMMG2D_int_nodes *)b;
  return (elementA->proc > elementB->proc) - (elementA->proc < elementB->proc);
}

int index_min(int *vec, int size) {
    int index_min = 0;

    for (int i = 1; i < size; i++) {
        if (vec[i] < vec[index_min]) {
            index_min = i;
        }
    }

    return index_min;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * Fill the list of interface nodes
 *
 */
void PMMG2D_fill_interface_nodes_list( PMMG2D_pParMesh parmesh )
{
  if (parmesh->nprocs == 1) return;

  int ier, size, k, n, size_local;
  int* list_local_interface_nodes = (int*) malloc(2*parmesh->mesh->np * sizeof(int));
  int rcounts[parmesh->nprocs];
  int displs[parmesh->nprocs];

  // Store the interface nodes locally with global node index
  size = 0;
  for ( k = 1; k <= parmesh->mesh->np; k++) {
    if ((&parmesh->mesh->point[k])->tag & MG_PARBDY) {
      list_local_interface_nodes[size] = k;
      list_local_interface_nodes[size+1] = (&parmesh->mesh->point[k])->ref;
      size += 2;
    }
  }

  size_local = size/2;

  PMMG2D_CALLOC(parmesh, parmesh->list_interface_nodes, size_local, PMMG2D_int_nodes, "interface nodes", ier = 0);
  InterfaceNode* list_interface_local_nodes = (InterfaceNode*) malloc(size_local * sizeof(InterfaceNode));
  for ( k = 0; k < size_local; k++) {
    list_interface_local_nodes[k].local_index = list_local_interface_nodes[2*k];
    list_interface_local_nodes[k].global_index = list_local_interface_nodes[2*k+1];
  }
  qsort(list_interface_local_nodes, size_local, sizeof(InterfaceNode), compare_interface_nodes);

  // Gather this list and broadcast it to every process
  MPI_CHECK( MPI_Allgather(&size, 1, MPI_INT, rcounts, 1, MPI_INT, parmesh->comm ), ier=6);

  displs[0] = 0;
  for ( k = 1; k < parmesh->nprocs; k++ ) {
    displs[k] = displs[k-1] + rcounts[k-1];
  }

  int* list_global_interface_nodes = (int*) malloc((displs[parmesh->nprocs-1]+rcounts[parmesh->nprocs-1]) * sizeof(int));
  MPI_CHECK( MPI_Allgatherv(&list_local_interface_nodes[0], size, MPI_INT,
                            list_global_interface_nodes, rcounts, displs,
                            MPI_INT, parmesh->comm), ier = 6);

  free(list_local_interface_nodes);

  if (size_local == 0) {
    free(list_global_interface_nodes);
    parmesh->nip = 0;
    return;
  }

  // Divide information in list_global_interface_nodes
  size = (displs[parmesh->nprocs-1]+rcounts[parmesh->nprocs-1])/2;
  InterfaceNode* list_interface_nodes = (InterfaceNode*) malloc(size * sizeof(InterfaceNode));
  int proc = 0;
  int i = 1;
  for ( k = 0; k < size; k++) {
    if (i < parmesh->nprocs) {
      if ( k == displs[i]/2 ) {
        i++;
        proc++;
        while ( i < parmesh->nprocs && rcounts[proc] == 0) {
          i++;
          proc++;
        }
      }
    }
    list_interface_nodes[k].proc = proc;
    list_interface_nodes[k].local_index = list_global_interface_nodes[2*k];
    list_interface_nodes[k].global_index = list_global_interface_nodes[2*k+1];
  }

  free(list_global_interface_nodes);

  // Sort the list relatively to the global_index
  qsort(list_interface_nodes, size, sizeof(InterfaceNode), compare_interface_nodes);

  // Fill the list of interface nodes
  i = 0;
  int l = 0;
  for ( k = 0; k < size; k++) {

    // Number of proc containing the global index
    n = 0;
    while (list_interface_nodes[k].global_index == list_interface_nodes[k+n].global_index) {
      n++;
      if (k+n == size) break;
    }

    // Check whether the node is within the processor domain
    int is_proc = 0;
    for (i = 0; i < n; i++) {
      if (list_interface_nodes[k+i].proc == parmesh->myrank) {
        is_proc = 1;
        break;
      }
    }

    if (!is_proc) continue;

    (&parmesh->list_interface_nodes[l])->index = (int*)malloc(n*sizeof(int));
    (&parmesh->list_interface_nodes[l])->proc = (int*)malloc(n*sizeof(int));
    (&parmesh->list_interface_nodes[l])->coord = (&parmesh->mesh->point[list_interface_local_nodes[i].local_index])->c;

    int index[n];
    int list_proc[n];
    int m = -1;
    for ( int j = 0; j < n; j++ ) {
      index[j] = list_interface_nodes[k+j].local_index;
      list_proc[j] = list_interface_nodes[k+j].proc;
      if (list_proc[j] == parmesh->myrank) m = j;
    }

    assert(m != -1);

    // The first component is the own process
    (&parmesh->list_interface_nodes[l])->index[0] = index[m];
    (&parmesh->list_interface_nodes[l])->proc[0] = parmesh->myrank;
    (&parmesh->list_interface_nodes[l])->nb = 1;
    list_proc[m] = parmesh->nprocs+1;

    // Then, fill in the ascending order relatively to the proc number
    for ( int j = 1; j < n; j++) {
      int m = index_min(list_proc,n);
      (&parmesh->list_interface_nodes[l])->index[j] = index[m];
      (&parmesh->list_interface_nodes[l])->proc[j] = list_proc[m];
      (&parmesh->list_interface_nodes[l])->nb++;
      list_proc[m] = parmesh->nprocs+100;
    }

    l++;
    k += n-1;
  }

  free(list_interface_nodes);
  free(list_interface_local_nodes);

  parmesh->nip = l;

}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * Free the list of interface nodes
 *
 */
void PMMG2D_free_interface_nodes_list ( PMMG2D_pParMesh parmesh )
{

  if (parmesh->list_interface_nodes != NULL) {
    for (int i = 0; i < parmesh->nip; ++i) {
      PMMG2D_pint_nodes node = &parmesh->list_interface_nodes[i];

      if (node != NULL) {

        if (node->index != NULL) {
          free(node->index);
        }
        if (node->proc != NULL) {
          free(node->proc);
        }
        node->nb = 0;
      }

    }

    parmesh->nip = 0;

    PMMG2D_DEL_MEM(parmesh, parmesh->list_interface_nodes, PMMG2D_int_nodes, "Delete list_interface_nodes");
  }

}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Update the mesh adjacency
 * Find the points on the parallel boundary between partitions and freeze them
 * Find the edges on the parallel boundary between partitions and freeze them
 *
 */
int PMMG2D_mark_parallel_interface_nodes( PMMG2D_pParMesh parmesh )
{
  int k, i, i1, i2;
  MMG5_pMesh   mesh;
  MMG5_pTria   pt;
  int     *adja;

  mesh = parmesh->mesh;

  // Compute the mesh adjacency
  if (mesh->adja) PMMG2D_DEL_MEM(mesh, mesh->adja, int, "remove adjacency table" );
  if ( (!mesh->adja) && mesh->nt && (MMG2D_hashTria( mesh ) != 1) ) {
    fprintf( stderr,"  ## PMMG2D Hashing problem (1).\n" );
    return 0;
  }

  for (k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) )  continue;
    adja = &mesh->adja[3*(k-1)+1];

    // Check whether there is an edge without neighbours
    for (i = 0; i < 3; i++) {

      if ( adja[i] ) continue; // there is a triangle neighbour

      if ( pt->tag[i] & MG_BDY ) continue; // it is a true boundary

      // Otherwise, the edge is a fictitious parallel boundary
      if ( (pt->tag[i] & MG_REQ) ) pt->tag[i] |= MG_PARBDYBDY;
      pt->tag[i] |= MG_PARBDY;

      // also mark the points
      i1  = MMG5_inxt2[i];
      i2  = MMG5_iprv2[i];
      (&mesh->point[pt->v[i1]])->tag |= MG_PARBDY;
      if ( (&mesh->point[pt->v[i1]])->tag & MG_REQ ) (&mesh->point[pt->v[i1]])->tag |= MG_PARBDYBDY;

      (&mesh->point[pt->v[i2]])->tag |= MG_PARBDY;
      if ( (&mesh->point[pt->v[i2]])->tag & MG_REQ ) (&mesh->point[pt->v[i2]])->tag |= MG_PARBDYBDY;
    }
  }

  for (k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];
    for (i = 0; i < 3; i++) {
      if (pt->tag[i] & MG_PARBDY) pt->tag[i] |= MG_REQ + MG_NOSURF;
    }
  }

  for (k = 1; k <= mesh->np; k++) {
    if ((&mesh->point[k])->tag & MG_PARBDY) (&mesh->point[k])->tag |= MG_REQ + MG_NOSURF;
  }

  return 1;
}

static double norm2( double x1, double y1, double x2, double y2 ) {
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}

// Write a vtk file containing the mesh and the metric to check
static void write_vtk_file(PMMG2D_pParMesh parmesh)
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
  int* is_boundary_nodes = (int*) malloc(parmesh->mesh->np * sizeof(int));
  for (int n = 0; n < parmesh->mesh->np; n++) is_boundary_nodes[n] = parmesh->mesh->point[n+1].tag;
  for (int n = 0; n < parmesh->nip; n++)
    is_boundary_nodes[(&parmesh->list_interface_nodes[n])->index[0]-1] = (&parmesh->list_interface_nodes[n])->proc[1]+1;

  // SCALARS fail when reading with paraview 5.10 for unknown reason
  fprintf(fp, "VECTORS Boundary int \n");
  for (int n = 0; n < parmesh->mesh->np; n++)
    fprintf(fp, "%d %d %d\n", is_boundary_nodes[n], 0, 0);

  fprintf(fp, "\n");

  free(is_boundary_nodes);

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
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Find the list of points at the interface between processors
 * Fill the list of interface nodes
 *
 */
int PMMG2D_bind_parallel_interfaces( PMMG2D_pParMesh parmesh )
{
  int k, j, n, i, npt_global, ier, idx, idy, closest_point_id;
  double x1, y1, x2, y2, w1, w2, h1, h2, min_dist;
  double* list_global_interface_nodes_x = NULL;
  double* list_global_interface_nodes_y = NULL;
  int* list_global_interface_nodes_id = NULL;
  double ESP_P_1 = 1. + PMMG2D_EPS;
  MMG5_pMesh   mesh;
  MPI_Request send_request, recv_request;
  MPI_Request send_request2, recv_request2;
  MPI_Status     status, status2;

  mesh = parmesh->mesh;

  // 1) Mark the parallel interface nodes as MG_PARBDY
  ier = PMMG2D_mark_parallel_interface_nodes( parmesh );

  // 2) Find the bounds of the interface nodes of each process
  n = 0;
  for (k = 1; k <= mesh->np; k++) {
    if ((&mesh->point[k])->tag & MG_PARBDY) n++;
  }

  double* list_interface_nodes = (double*) malloc(2*n*sizeof(double));
  int* list_interface_nodes_id = (int*) malloc(n*sizeof(int));
  n = 0;
  for (k = 1; k <= mesh->np; k++) {
    if ((&mesh->point[k])->tag & MG_PARBDY) {
      list_interface_nodes[2*n] = (&mesh->point[k])->c[0];
      list_interface_nodes[2*n+1] = (&mesh->point[k])->c[1];
      list_interface_nodes_id[n] = k;
      n++;
    }
  }

  double minX = list_interface_nodes[0];
  double maxX = list_interface_nodes[0];
  double minY = list_interface_nodes[1];
  double maxY = list_interface_nodes[1];

  for (k = 1; k < n; k++) {
    if (list_interface_nodes[2*k] < minX) minX = list_interface_nodes[2*k];
    if (list_interface_nodes[2*k] > maxX) maxX = list_interface_nodes[2*k];
    if (list_interface_nodes[2*k+1] < minY) minY = list_interface_nodes[2*k+1];
    if (list_interface_nodes[2*k+1] > maxY) maxY = list_interface_nodes[2*k+1];
  }

  // 3) Gather the bounds over the processors and associate the neighbours
  double list_minX[parmesh->nprocs], list_maxX[parmesh->nprocs];
  double list_minY[parmesh->nprocs], list_maxY[parmesh->nprocs];
  MPI_CHECK( MPI_Allgather(&minX, 1, MPI_DOUBLE, list_minX, 1, MPI_DOUBLE, parmesh->comm), ier = 0 );
  MPI_CHECK( MPI_Allgather(&minY, 1, MPI_DOUBLE, list_minY, 1, MPI_DOUBLE, parmesh->comm), ier = 0 );
  MPI_CHECK( MPI_Allgather(&maxX, 1, MPI_DOUBLE, list_maxX, 1, MPI_DOUBLE, parmesh->comm), ier = 0 );
  MPI_CHECK( MPI_Allgather(&maxY, 1, MPI_DOUBLE, list_maxY, 1, MPI_DOUBLE, parmesh->comm), ier = 0 );

  int list_neighbours_proc[parmesh->nprocs];
  int number_points[parmesh->nprocs];
  for (k = 0; k < parmesh->nprocs; k++) {
    list_neighbours_proc[k] = -1;
    number_points[k] = 0;
  }
    
  x1 = minX;
  y1 = minY;
  w1 = maxX-minX;
  h1 = maxY-minY;

  i = 0;
  for (k = 0; k < parmesh->nprocs; k++) {

    if (k == parmesh->myrank) continue;

    x2 = list_minX[k];
    y2 = list_minY[k];
    w2 = list_maxX[k] - list_minX[k];
    h2 = list_maxY[k] - list_minY[k];

    // Check whether the rectangles intersect
    if ( (x1 <= x2+ESP_P_1*w2) && (x2 <= x1+ESP_P_1*w1) && (y1 <= y2+ESP_P_1*h2) && (y2 <= y1+ESP_P_1*h1) ) {
      list_neighbours_proc[i] = k;
      i++;
    }

  }

  // 4) Exchange the interface nodes between the neighbouring processors
  i = 0;
  npt_global = 0;
  while (list_neighbours_proc[i] != -1) {

    k = list_neighbours_proc[i];

    // Send the number of interface nodes
    MPI_CHECK( MPI_Send(&n, 1, MPI_INT, k, (parmesh->myrank+1)*MPI_ANALYS_TAG + k, parmesh->comm), ier = 0 );

    // Receive the number of interface nodes
    MPI_CHECK( MPI_Recv(&number_points[k], 1, MPI_INT, k, (k+1)*MPI_ANALYS_TAG + parmesh->myrank,
                         parmesh->comm, &status), ier = 0 );

    // Send the list of interface nodes
    MPI_CHECK( MPI_Isend(&list_interface_nodes[0], 2*n, MPI_DOUBLE, k,
                         (parmesh->myrank+1)*MPI_ANALYS_TAG + parmesh->nprocs + 1 + k, 
                         parmesh->comm, &send_request), ier = 0 );
    MPI_CHECK( MPI_Isend(&list_interface_nodes_id[0], n, MPI_INT, k,
                         (parmesh->myrank+1)*MPI_ANALYS_TAG + 2*parmesh->nprocs + 1 + k,
                         parmesh->comm, &send_request2), ier = 0 );

    npt_global += number_points[k];
    double* local_list = (double*) malloc(2*number_points[k] * sizeof(double));
    int* local_list_id = (int*) malloc(number_points[k] * sizeof(int));

    MPI_CHECK( MPI_Irecv(local_list, 2*number_points[k], MPI_DOUBLE, k,
                         (k+1)*MPI_ANALYS_TAG + parmesh->myrank + parmesh->nprocs + 1,
                         parmesh->comm, &recv_request), ier = 0 );
    MPI_CHECK( MPI_Irecv(local_list_id, number_points[k], MPI_INT, k,
                         (k+1)*MPI_ANALYS_TAG + parmesh->myrank + 2*parmesh->nprocs + 1,
                         parmesh->comm, &recv_request2), ier = 0 );

    MPI_Wait(&send_request, &status);
    MPI_Wait(&recv_request, &status);
    MPI_Wait(&send_request2, &status2);
    MPI_Wait(&recv_request2, &status2);

    double* local_list_x = (double*) malloc(number_points[k] * sizeof(double));
    double* local_list_y = (double*) malloc(number_points[k] * sizeof(double));

    for (int l = 0; l < number_points[k]; l++) {
      local_list_x[l] = local_list[2*l];
      local_list_y[l] = local_list[2*l+1];
    }

    PMMG2D_REALLOC(mesh, list_global_interface_nodes_x, npt_global, npt_global-number_points[k],
                   double, "Exchange interface nodes in PMMG2D_bind_parallel_interfaces",return 0);
    memcpy(list_global_interface_nodes_x + npt_global - number_points[k], 
           local_list_x, number_points[k]*sizeof(double));

    PMMG2D_REALLOC(mesh, list_global_interface_nodes_y, npt_global, npt_global-number_points[k],
                   double, "Exchange interface nodes in PMMG2D_bind_parallel_interfaces",return 0);
    memcpy(list_global_interface_nodes_y + npt_global - number_points[k], 
           local_list_y, number_points[k]*sizeof(double));

    PMMG2D_REALLOC(mesh, list_global_interface_nodes_id, npt_global, npt_global-number_points[k],
                   int, "Exchange interface nodes in PMMG2D_bind_parallel_interfaces",return 0);
    memcpy(list_global_interface_nodes_id + npt_global - number_points[k], 
           local_list_id, number_points[k]*sizeof(int));

    i++;

    free(local_list);
    free(local_list_id);
    free(local_list_x);
    free(local_list_y);
  }

  // 5) Use k-d tree algorithm to find the closest points 
  // First build a grid and store the neighbour points in each cell
  int nt_tot;
  MPI_Allreduce(&mesh->nt, &nt_tot, 1, MPI_INT, MPI_SUM, parmesh->comm);
  int GRID_SIZE = sqrt(nt_tot/50);
  double cellWidth = (maxX - minX) / GRID_SIZE;
  double cellHeight = (maxY - minY) / GRID_SIZE;

  int*** grid = (int***) malloc(GRID_SIZE * sizeof(int**));
  for (i = 0; i < GRID_SIZE; i++) {
    grid[i] = (int**) malloc(GRID_SIZE * sizeof(int*));
    for (k = 0; k < GRID_SIZE; k++) grid[i][k] = (int*) malloc(2*sizeof(int));
  }
  int** counts = (int**) malloc(GRID_SIZE * sizeof(int*));
  for (i = 0; i < GRID_SIZE; i++) counts[i] = (int*) malloc(GRID_SIZE * sizeof(int));

  int list_idx[2], list_idy[2];
  for (i = 0; i < GRID_SIZE; i++) {for (k = 0; k < GRID_SIZE; k++) counts[i][k] = 0;}

  j = number_points[list_neighbours_proc[0]];
  k = 0;
  for (i = 0; i < npt_global; i++) {

    if (i == j) {
      k++;
      j += number_points[list_neighbours_proc[k]];
    }

    // For the sake of security, if a point is on a grid boundary, 
    // it must belong to all the surrounding cells
    list_idx[0] = (int)((list_global_interface_nodes_x[i] - minX) / cellWidth + PMMG2D_EPS);
    list_idx[1] = (int)((list_global_interface_nodes_x[i] - minX) / cellWidth - PMMG2D_EPS);
    if (list_idx[1] == list_idx[0]) list_idx[1] = -1;
    list_idy[0] = (int)((list_global_interface_nodes_y[i] - minY) / cellHeight + PMMG2D_EPS);
    list_idy[1] = (int)((list_global_interface_nodes_y[i] - minY) / cellHeight - PMMG2D_EPS);
    if (list_idy[1] == list_idy[0]) list_idy[1] = -1;

    for (int lx = 0; lx < 2; lx++) {
      for (int ly = 0; ly < 2; ly++) {
        idx = list_idx[lx];
        idy = list_idy[ly];
        // Check whether the vertex is inside the grid domain
        if (idx >= 0 && idx < GRID_SIZE && idy >= 0 && idy < GRID_SIZE) {
          grid[idx][idy][2*counts[idx][idy]] = list_neighbours_proc[k];
          grid[idx][idy][2*counts[idx][idy]+1] = i;
          counts[idx][idy]++;
          grid[idx][idy] = (int*)realloc(grid[idx][idy], 2*(counts[idx][idy]+1)*sizeof(int));
        }
      }
    }
  }

  // Second associate the local points to the neighbour points in grid
  PMMG2D_pint_nodes interface_nodes;
  PMMG2D_CALLOC(parmesh, interface_nodes, n, PMMG2D_int_nodes, "interface nodes (2)", ier = 0);

  for (i = 0; i < n; i++) {
    (&interface_nodes[i])->index = (int*)malloc(sizeof(int));
    (&interface_nodes[i])->proc = (int*)malloc(sizeof(int));
    (&interface_nodes[i])->coord = (double*) malloc(2*sizeof(double));
    (&interface_nodes[i])->coord[0] = list_interface_nodes[2*i];
    (&interface_nodes[i])->coord[1] = list_interface_nodes[2*i+1];

    // The first component is the own process
    (&interface_nodes[i])->index[0] = list_interface_nodes_id[i];
    (&interface_nodes[i])->proc[0] = parmesh->myrank;
    (&interface_nodes[i])->nb = 1;
  }

  for (i = 0; i < n; i++) {
    idx = (int)((list_interface_nodes[2*i] - minX) / cellWidth);
    // If the point is on a boundary, idx can be outside the domain (<0 or >= GRID_SIZE)
    // It must be corrected
    if (idx == -1) idx = 0;
    if (idx == GRID_SIZE) idx = GRID_SIZE-1;
    idy = (int)((list_interface_nodes[2*i+1] - minY) / cellHeight);
    if (idy == -1) idy = 0;
    if (idy == GRID_SIZE) idy = GRID_SIZE-1;

    for (k = 0; k < counts[idx][idy]; k++) {
      if (norm2(list_interface_nodes[2*i], list_interface_nodes[2*i+1], 
                list_global_interface_nodes_x[grid[idx][idy][2*k+1]], list_global_interface_nodes_y[grid[idx][idy][2*k+1]]) < PMMG2D_EPS*cellHeight*cellWidth) {
        (&interface_nodes[i])->index[(&interface_nodes[i])->nb] = list_global_interface_nodes_id[grid[idx][idy][2*k+1]];
        (&interface_nodes[i])->proc[(&interface_nodes[i])->nb] = grid[idx][idy][2*k];
        (&interface_nodes[i])->nb++;
      }
    }

    assert( (&interface_nodes[i])->nb != 1 );

  }

  free(list_interface_nodes);
  free(list_interface_nodes_id);

  // 6) Renumber the nodes with the correct global reference ids
  PMMG2D_renumber_distributed_mesh( parmesh, interface_nodes, n);

  // 7) Update the list of interface nodes 
  PMMG2D_fill_interface_nodes_list( parmesh );

  // To check
  //write_vtk_file(parmesh);

  for (i = 0; i < GRID_SIZE; i++) {
    for (k = 0; k < GRID_SIZE; k++) free(grid[i][k]);
    free(grid[i]);
    free(counts[i]);
  }
  free(grid);
  free(counts);

  PMMG2D_DEL_MEM(mesh, list_global_interface_nodes_x, double, "Delete interface nodes in PMMG2D_bind_parallel_interfaces");
  PMMG2D_DEL_MEM(mesh, list_global_interface_nodes_y, double, "Delete interface nodes in PMMG2D_bind_parallel_interfaces");
  PMMG2D_DEL_MEM(mesh, list_global_interface_nodes_id, int, "Delete interface nodes in PMMG2D_bind_parallel_interfaces");

  for (i = 0; i < n; i++) {
    if ( NULL != (&interface_nodes[i])->index ) {
      free((&interface_nodes[i])->index);
    }
    if ( NULL != (&interface_nodes[i])->proc ) {
      free((&interface_nodes[i])->proc);
    }
    (&interface_nodes[i])->nb = 0;
  }

  PMMG2D_DEL_MEM(parmesh, interface_nodes, PMMG2D_int_nodes, "Delete interface nodes");

  return 1;
}

static int min(int* list) {
  int res = list[0];
  if (list[1] < res) res = list[1];
  return res;
}

static int max(int* list) {
  int res = list[0];
  if (list[1] > res) res = list[1];
  return res;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Correct the MG_BDY tag at the interface between processors
 * That could have been removed in PMMG2D_preprocessMesh function
 *
 */
int PMMG2D_correct_BDY_tags( PMMG2D_pParMesh parmesh )
{
  int k,i;
  MMG5_pTria pt;
  MMG5_pMesh mesh;

  mesh = parmesh->mesh;

  int ier = MMG5_unscaleMesh(mesh,parmesh->met,NULL) ;

  // Correct the NOSURF tags removed by MMG2D_analysis improperly
  for (int k = 1; k <= parmesh->mesh->np; k++ ) {
    if ( parmesh->mesh->info.nosurf && (parmesh->mesh->point[k].tag & MG_BDY) && !(parmesh->mesh->point[k].tag & MG_REQ) ) {
      parmesh->mesh->point[k].tag |= MG_REQ + MG_NOSURF;
    }
  }
 
  // Correct the MG_BDY tags imposed on parallel interface improperly in MMG2D_analysis
  int* n_points_bdy = (int*) malloc(mesh->np * sizeof(int));
  int* list_remove_points = (int*) malloc(mesh->np * sizeof(int));
  int* list_remove_points_no_surf = (int*) malloc(mesh->np * sizeof(int));
  int n_remove_points = 0;
  for (k = 0; k < mesh->np; k++) {
    n_points_bdy[k] = 0;
    list_remove_points_no_surf[k] = 0;
  }

  for (k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {

      int p1 = pt->v[MMG5_inxt2[i]];
      int p2 = pt->v[MMG5_iprv2[i]];
      if (pt->edg[i] == -1) {
        if (pt->tag[i] & MG_BDY) {
          n_points_bdy[pt->v[MMG5_inxt2[i]]-1]++;
          n_points_bdy[pt->v[MMG5_iprv2[i]]-1]++;
        }

        pt->tag[i] &= ~MG_BDY;
        pt->tag[i] &= ~MG_CRN;
        pt->tag[i] &= ~MG_REF;
        pt->tag[i] &= ~MG_GEO;
        if (pt->tag[i] & MG_NOSURF) {
          pt->tag[i] &= ~MG_REQ;
          pt->tag[i] &= ~MG_NOSURF;
        }

        if (n_points_bdy[pt->v[MMG5_inxt2[i]]-1] == 2) {
          mesh->point[pt->v[MMG5_inxt2[i]]].tag &= ~MG_BDY;
          mesh->point[pt->v[MMG5_inxt2[i]]].tag &= ~MG_REF;
          mesh->point[pt->v[MMG5_inxt2[i]]].tag &= ~MG_GEO;
          mesh->point[pt->v[MMG5_inxt2[i]]].tag &= ~MG_CRN;
          if (mesh->point[pt->v[MMG5_inxt2[i]]].tag & MG_NOSURF) {
            mesh->point[pt->v[MMG5_inxt2[i]]].tag &= ~MG_NOSURF;
            mesh->point[pt->v[MMG5_inxt2[i]]].tag &= ~MG_REQ;
            list_remove_points_no_surf[n_remove_points] = 1;
          }
          list_remove_points[n_remove_points++] = pt->v[MMG5_inxt2[i]];
        }
        if (n_points_bdy[pt->v[MMG5_iprv2[i]]-1] == 2) {
          mesh->point[pt->v[MMG5_iprv2[i]]].tag &= ~MG_BDY;
          mesh->point[pt->v[MMG5_iprv2[i]]].tag &= ~MG_REF;
          mesh->point[pt->v[MMG5_iprv2[i]]].tag &= ~MG_GEO;
          mesh->point[pt->v[MMG5_iprv2[i]]].tag &= ~MG_CRN;
          if (mesh->point[pt->v[MMG5_iprv2[i]]].tag & MG_NOSURF) {
            mesh->point[pt->v[MMG5_iprv2[i]]].tag &= ~MG_NOSURF;
            mesh->point[pt->v[MMG5_iprv2[i]]].tag &= ~MG_REQ;
            list_remove_points_no_surf[n_remove_points] = 1;
          }
          list_remove_points[n_remove_points++] = pt->v[MMG5_iprv2[i]];
        }
      }
    }
  }

  // The case of an isolated MG_BDY point in local meshes needs further treatment
  int n = 0;
  for (k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) {

      if ( !(pt->tag[i] & MG_BDY) ) continue;

      if (!(mesh->point[pt->v[MMG5_inxt2[i]]].tag & MG_BDY)) {

        while (list_remove_points[n] != pt->v[MMG5_inxt2[i]]) n++;

        mesh->point[pt->v[MMG5_inxt2[i]]].tag |= MG_BDY;
        mesh->point[pt->v[MMG5_inxt2[i]]].tag |= MG_REF;
        mesh->point[pt->v[MMG5_inxt2[i]]].tag |= MG_GEO;
        mesh->point[pt->v[MMG5_inxt2[i]]].tag |= MG_CRN;
        if (list_remove_points_no_surf[n]) {
          mesh->point[pt->v[MMG5_inxt2[i]]].tag |= MG_NOSURF;
          mesh->point[pt->v[MMG5_inxt2[i]]].tag |= MG_REQ;
        }

      }

      if (!(mesh->point[pt->v[MMG5_iprv2[i]]].tag & MG_BDY)) {

        while (list_remove_points[n] != pt->v[MMG5_iprv2[i]]) n++;

        mesh->point[pt->v[MMG5_iprv2[i]]].tag |= MG_BDY;
        mesh->point[pt->v[MMG5_iprv2[i]]].tag |= MG_REF;
        mesh->point[pt->v[MMG5_iprv2[i]]].tag |= MG_GEO;
        mesh->point[pt->v[MMG5_iprv2[i]]].tag |= MG_CRN;
        if (list_remove_points_no_surf[n]) {
          mesh->point[pt->v[MMG5_iprv2[i]]].tag |= MG_NOSURF;
          mesh->point[pt->v[MMG5_iprv2[i]]].tag |= MG_REQ;
        }

      }

    }

  }

  // Find the bounding box
  double x_min = 1.e20, x_max = -1.e20, y_min = 1.e20, y_max = -1.e20;
  double* list_bdy = (double*) malloc(2*mesh->np * sizeof(double));
  n = 0;
  for (k = 1; k <= mesh->np; k++) {
    if (x_min > mesh->point[k].c[0]) x_min = mesh->point[k].c[0];
    if (y_min > mesh->point[k].c[1]) y_min = mesh->point[k].c[1];
    if (x_max < mesh->point[k].c[0]) x_max = mesh->point[k].c[0];
    if (y_max < mesh->point[k].c[1]) y_max = mesh->point[k].c[1];
    if (mesh->point[k].tag & MG_BDY) {
      list_bdy[n++] = mesh->point[k].c[0];
      list_bdy[n++] = mesh->point[k].c[1];
    }
  }

  int size;
  MPI_Allreduce(&n, &size, 1, MPI_INT, MPI_SUM, parmesh->comm);

  int rcounts[parmesh->nprocs];
  int displs[parmesh->nprocs];
  MPI_Allgather(&n, 1, MPI_INT, rcounts, 1, MPI_INT, parmesh->comm);

  displs[0] = 0;
  for (i = 1; i < parmesh->nprocs; i++) displs[i] = displs[i-1] + rcounts[i-1];

  double* list_bdy_global = (double*) malloc(size * sizeof(double));
  MPI_Allgatherv(&list_bdy[0], n, MPI_DOUBLE, list_bdy_global, rcounts, displs, MPI_DOUBLE, parmesh->comm);

  // Use k-d tree algorithm to find the closest points 
  // First build a grid and store the boundary points in each cell
  double maxX, maxY, minX, minY;
  MPI_Allreduce(&x_min, &minX, 1, MPI_DOUBLE, MPI_MIN, parmesh->comm);
  MPI_Allreduce(&y_min, &minY, 1, MPI_DOUBLE, MPI_MIN, parmesh->comm);
  MPI_Allreduce(&x_max, &maxX, 1, MPI_DOUBLE, MPI_MAX, parmesh->comm);
  MPI_Allreduce(&y_max, &maxY, 1, MPI_DOUBLE, MPI_MAX, parmesh->comm);

  if (!n_remove_points) {
    free(list_remove_points_no_surf);
    free(list_remove_points);
    free(n_points_bdy);
    free(list_bdy);
    free(list_bdy_global);
    return 1;
  }

  int nt_tot;
  MPI_Allreduce(&mesh->nt, &nt_tot, 1, MPI_INT, MPI_SUM, parmesh->comm);
  int GRID_SIZE = sqrt(nt_tot/50);
  double cellWidth = (maxX - minX) / GRID_SIZE;
  double cellHeight = (maxY - minY) / GRID_SIZE;

  int*** grid = (int***) malloc(GRID_SIZE * sizeof(int**)); 
  for (i = 0; i < GRID_SIZE; i++) {
    grid[i] = (int**) malloc(GRID_SIZE * sizeof(int*));
    for (k = 0; k < GRID_SIZE; k++) grid[i][k] = (int*) malloc(sizeof(int));
  }

  int** counts = (int**) malloc(GRID_SIZE * sizeof(int*));
  for (i = 0; i < GRID_SIZE; i++) counts[i] = (int*) malloc(GRID_SIZE * sizeof(int));
  int list_idx[2], list_idy[2];
  for (i = 0; i < GRID_SIZE; i++) {for (k = 0; k < GRID_SIZE; k++) counts[i][k] = 0;}

  for (i = 0; i < size/2; i++) {
    // For the sake of security, if a point is on a grid boundary, 
    // it must belong to all the surrounding cells
    list_idx[0] = (int)((list_bdy_global[2*i] - minX) / cellWidth + PMMG2D_EPS);
    list_idx[1] = (int)((list_bdy_global[2*i] - minX) / cellWidth - PMMG2D_EPS);
    int minx = min(list_idx);
    if (minx < 0) minx = 0;
    int maxx = max(list_idx);
    if (maxx >= GRID_SIZE) maxx = GRID_SIZE-1;
    list_idy[0] = (int)((list_bdy_global[2*i+1] - minY) / cellHeight + PMMG2D_EPS);
    list_idy[1] = (int)((list_bdy_global[2*i+1] - minY) / cellHeight - PMMG2D_EPS);
    int miny = min(list_idy);
    if (miny < 0) miny = 0;
    int maxy = max(list_idy);
    if (maxy >= GRID_SIZE) maxy = GRID_SIZE-1;

    for (int idx = minx; idx <= maxx; idx++) {
      for (int idy = miny; idy <= maxy; idy++) {
        grid[idx][idy][counts[idx][idy]] = 2*i;
        counts[idx][idy]++;
        grid[idx][idy] = (int*)realloc(grid[idx][idy], (counts[idx][idy]+1)*sizeof(int));
      }
    }
  }

  double element_size = sqrt(pow(mesh->point[mesh->tria[1].v[MMG5_inxt2[0]]].c[0] - mesh->point[mesh->tria[1].v[MMG5_iprv2[0]]].c[0],2)
                            +pow(mesh->point[mesh->tria[1].v[MMG5_inxt2[0]]].c[1] - mesh->point[mesh->tria[1].v[MMG5_iprv2[0]]].c[1],2));

  for (i = 0; i < n_remove_points; i++) {

    int id = list_remove_points[i];

    list_idx[0] = (int)((mesh->point[id].c[0] - minX) / cellWidth + PMMG2D_EPS);
    list_idx[1] = (int)((mesh->point[id].c[0] - minX) / cellWidth - PMMG2D_EPS);
    int minx = min(list_idx);
    if (minx < 0) minx = 0;
    int maxx = max(list_idx);
    if (maxx >= GRID_SIZE) maxx = GRID_SIZE-1;

    list_idy[0] = (int)((mesh->point[id].c[1] - minY) / cellHeight + PMMG2D_EPS);
    list_idy[1] = (int)((mesh->point[id].c[1] - minY) / cellHeight - PMMG2D_EPS);
    int miny = min(list_idy);
    if (miny < 0) miny = 0;
    int maxy = max(list_idy);
    if (maxy >= GRID_SIZE) maxy = GRID_SIZE-1;

    for (int idx = minx; idx <= maxx; idx++) {
      for (int idy = miny; idy <= maxy; idy++) {
        for (n = 0; n < counts[idx][idy]; n++) {

          if (fabs(mesh->point[id].c[0] - list_bdy_global[grid[idx][idy][n]]) < element_size*PMMG2D_EPS &&
              fabs(mesh->point[id].c[1] - list_bdy_global[grid[idx][idy][n]+1]) < element_size*PMMG2D_EPS) {

            // Restore the tags that have been removed incorrectly
            mesh->point[id].tag |= MG_BDY;
            mesh->point[id].tag |= MG_REF;
            mesh->point[id].tag |= MG_GEO;
            mesh->point[id].tag |= MG_CRN;
            if (list_remove_points_no_surf[i]) {
              mesh->point[id].tag |= MG_NOSURF;
              mesh->point[id].tag |= MG_REQ;
            }

            break;

          }

        }

      }
    }
  }

  free(list_bdy_global);
  free(list_remove_points_no_surf);
  free(list_remove_points);
  free(n_points_bdy);
  free(list_bdy);
  for (i = 0; i < GRID_SIZE; i++) {
    for (k = 0; k < GRID_SIZE; k++) free(grid[i][k]);
    free(grid[i]);
    free(counts[i]);
  }
  free(grid);
  free(counts);

  ier = MMG5_scaleMesh(mesh,parmesh->met,NULL);

  return ier;
}
