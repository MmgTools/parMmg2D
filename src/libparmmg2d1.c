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
 * \file libparmmg2d1.c
 * \brief Wrapper for the parallel remeshing library.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Internal function that perform the parallel remeshing.
 *
 */
#include "parmmg2d.h"

// Write a vtk file containing the mesh and the metric to check
static void check_vtk_file(PMMG2D_pParMesh parmesh)
{
  MMG5_pTria pt;
  char name_file[20];
  sprintf(name_file, "remesh_%d.vtk", parmesh->myrank);

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
  for (int n = 0; n < parmesh->mesh->np; n++) {
    if (parmesh->mesh->point[n+1].tag & MG_PARBDY) {
        is_boundary_nodes[n] = -1;//parmesh->mesh->point[n].ref;
    }
    else {
        is_boundary_nodes[n] = 0;
    }
  }

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
      fprintf(fp, "%.14f %.14f %.14f\n", parmesh->met->m[3*n], parmesh->met->m[3*n+1], parmesh->met->m[3*n+2]);
    }
  }

  fclose(fp);
}

/**
 * \param parmesh pointer toward a parmesh structure where the boundary entities
 * are stored into xtetra and xpoint strucutres
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the initial meshes and metrics
 *
**/
int PMMG2D_copy_initial_mesh_met( PMMG2D_pParMesh parmesh )
{
  MMG5_pTria      pt,ptCur;
  MMG5_pPoint     ppt,pptCur;
  int             *adja,*oldAdja;
  int i;

  parmesh->initial_mesh = NULL;
  parmesh->initial_met  = NULL;

  // If there is not mesh, nothing to copy
  if (!parmesh->mesh->nt) return 1;

  MMG2D_Init_mesh( MMG5_ARG_start,
                   MMG5_ARG_ppMesh, &parmesh->initial_mesh,
                   MMG5_ARG_ppMet, &parmesh->initial_met,
                   MMG5_ARG_end );

  // Set maximum memory
  parmesh->initial_mesh->memMax = parmesh->memGloMax;

  // Set sizes
  if ( !PMMG2D_setMeshSize( parmesh->initial_mesh, parmesh->mesh->np, parmesh->mesh->nt,0,0 ) ) return 0;

  PMMG2D_CALLOC(parmesh->initial_mesh, parmesh->initial_mesh->adja, 3*parmesh->initial_mesh->ntmax+4,int,"tria adjacency table",return 0);

  // Set metric size
  if ( parmesh->info.inputMet == 1 ) {
    if ( !MMG2D_Set_solSize(parmesh->initial_mesh,parmesh->initial_met,MMG5_Vertex,parmesh->mesh->np,parmesh->met->type) )
      return 0;
  }

  // Copy the info structure of the initial mesh: it contains the remeshing options
  if ( !PMMG2D_copy_mmgInfo ( &parmesh->mesh->info, &parmesh->initial_mesh->info ) ) return 0;

  // Loop on trias
  for ( i = 1; i < parmesh->mesh->nt+1; ++i ) {
    pt = &parmesh->mesh->tria[i];
    ptCur = &parmesh->initial_mesh->tria[i];

    if ( !MG_EOK(pt) ) continue;

    // Copy tria
    memcpy( ptCur, pt, sizeof(MMG5_Tria) );

    // Copy element's adjacency
    assert( parmesh->mesh->adja );
    if( parmesh->mesh->adja ) {
      oldAdja    =    &parmesh->initial_mesh->adja[ 3*( i-1 )+1 ];
      adja = &parmesh->mesh->adja[ 3*( i-1 )+1 ];
      memcpy( oldAdja, adja, 3*sizeof(int) );
    }

  }

  // Loop on points
  for ( i = 1; i < parmesh->mesh->np+1; ++i ) {
    ppt = &parmesh->mesh->point[i];
    pptCur = &parmesh->initial_mesh->point[i];

    if ( !MG_VOK(ppt) ) {

      // Only copy the tag (to detect the not VOK point)
      pptCur->tag = ppt->tag;

    } else {

      // Copy point
      memcpy( pptCur, ppt, sizeof(MMG5_Point) );

      // Copy metrics
      if ( parmesh->info.inputMet == 1 ) {
        memcpy( &parmesh->initial_met->m[ i*parmesh->met->size ], &parmesh->met->m[i*parmesh->met->size], parmesh->met->size*sizeof(double) );
      }

    }
  }

  return 1;

}

/**
 * \param parmesh pointer toward a parmesh structure where the boundary entities
 * are stored into xtetra and xpoint strucutres
 *
 * \return 0 if fail, 1 if success
 *
 * Copy the old meshes and metrics
 *
**/
int PMMG2D_copy_mesh_met( PMMG2D_pParMesh parmesh )
{
  MMG5_pTria      pt,ptCur;
  MMG5_pPoint     ppt,pptCur;
  int             *adja,*oldAdja;
  int i;

  parmesh->old_mesh = NULL;
  parmesh->old_met  = NULL;

  // If there is not mesh, nothing to copy
  if (!parmesh->mesh->nt) return 1;

  MMG2D_Init_mesh( MMG5_ARG_start,
                   MMG5_ARG_ppMesh, &parmesh->old_mesh,
                   MMG5_ARG_ppMet, &parmesh->old_met,
                   MMG5_ARG_end );

  // Set maximum memory
  parmesh->old_mesh->memMax = parmesh->memGloMax;

  // Set sizes
  if ( !PMMG2D_setMeshSize( parmesh->old_mesh, parmesh->mesh->np, parmesh->mesh->nt,0,0 ) ) return 0;

  PMMG2D_CALLOC(parmesh->old_mesh, parmesh->old_mesh->adja, 3*parmesh->old_mesh->ntmax+4,int,"tria adjacency table",return 0);

  // Set metric size
  if ( parmesh->info.inputMet == 1 ) {
    if ( !MMG2D_Set_solSize(parmesh->old_mesh,parmesh->old_met,MMG5_Vertex,parmesh->mesh->np,parmesh->met->type) )
      return 0;
  }

  // Copy the info structure of the initial mesh: it contains the remeshing options
  if ( !PMMG2D_copy_mmgInfo ( &parmesh->mesh->info, &parmesh->old_mesh->info ) ) return 0;

  // Loop on trias
  for ( i = 1; i < parmesh->mesh->nt+1; ++i ) {
    pt = &parmesh->mesh->tria[i];
    ptCur = &parmesh->old_mesh->tria[i];

    if ( !MG_EOK(pt) ) continue;

    // Copy tria
    memcpy( ptCur, pt, sizeof(MMG5_Tria) );

    // Copy element's adjacency
    assert( parmesh->mesh->adja );
    if( parmesh->mesh->adja ) {
      oldAdja    =    &parmesh->old_mesh->adja[ 3*( i-1 )+1 ];
      adja = &parmesh->mesh->adja[ 3*( i-1 )+1 ];
      memcpy( oldAdja, adja, 3*sizeof(int) );
    }

  }

  // Loop on points
  for ( i = 1; i < parmesh->mesh->np+1; ++i ) {
    ppt = &parmesh->mesh->point[i];
    pptCur = &parmesh->old_mesh->point[i];

    if ( !MG_VOK(ppt) ) {

      // Only copy the tag (to detect the not VOK point)
      pptCur->tag = ppt->tag;

    } else {

      // Copy point
      memcpy( pptCur, ppt, sizeof(MMG5_Point) );

      // Copy metrics
      if ( parmesh->info.inputMet == 1 ) {
        memcpy( &parmesh->old_met->m[ i*parmesh->met->size ], &parmesh->met->m[i*parmesh->met->size], parmesh->met->size*sizeof(double) );
      }

    }
  }

  return 1;

}

/**
 * \param parmesh pointer toward a parmesh structure where the boundary entities
 * are stored into xtetra and xpoint strucutres
 *
 * Main program of the parallel remeshing library: split the meshes over each
 * proc into groups, then perform niter of sequential remeshing of each group
 * (with moving of the proc boundaries between two iterations) and last, merge
 * the groups over each proc.
 *
 * \return PMMG2D_STRONGFAILURE if  we can't save the mesh (non-conform),
 *         PMMG2D_LOWFAILURE    if  we can save the mesh
 *         PMMG2D_SUCCESS
 */
int PMMG2D_parmmg2dlib1( PMMG2D_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  MMG5_pTria pt;
  mytime     ctim[TIMEMAX];
  int        ier,ier_end,ieresult,k,np_tot,np;
  int8_t     tim;
  char       stim[32];

  tminit(ctim,TIMEMAX);

  ier_end = PMMG2D_SUCCESS;

  assert ( parmesh->mesh );

  /** Set inputMet flag */
  parmesh->info.inputMet = 0;
  if ( parmesh->met && parmesh->met->m ) parmesh->info.inputMet = 1;

  ier = 1;
#ifndef NDEBUG
  inputMet = 0;
  MPI_CHECK( MPI_Allreduce( &parmesh->info.inputMet,&inputMet,1,MPI_UNSIGNED_CHAR,MPI_MAX,
                            parmesh->comm ),ier = 0 );

  if ( inputMet != parmesh->info.inputMet ) {
    printf ("  ## Warning: input metric not provided on rank %d while provided on others.\n", parmesh->myrank);
    parmesh->info.inputMet = inputMet;
  }
#endif

  // If asked, compute the contour of the partitions and copy the initial mesh
  int size;
  double **polygon = NULL;

  if (parmesh->info.optim_interp) {
    polygon = build_partition_contour(parmesh->mesh, &size);
    if ( !PMMG2D_copy_initial_mesh_met( parmesh )) {
      fprintf(stderr,"\n  ## Error: %s: unable to copy initial mesh\n",__func__);
      return 0;
    }
  }

  // Mesh adaptation loop
  for ( parmesh->iter = 0; parmesh->iter < parmesh->niter; parmesh->iter++ ) {

    if ( parmesh->info.imprim > PMMG2D_VERB_STEPS ) {
      tim = 1;
      if ( parmesh->iter > 0 ) {
        chrono(OFF,&(ctim[tim]));
      }
      if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
        fprintf(stdout,"\n" );
      }

      printim(ctim[tim].gdif,stim);
      chrono(ON,&(ctim[tim]));
      fprintf(stdout,"\r       adaptation: iter %d   cumul. timer %s",parmesh->iter+1,stim);fflush(stdout);
    }

    tim = 4;
    if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
      chrono(RESET,&(ctim[tim]));
      chrono(ON,&(ctim[tim]));
    }

    mesh = parmesh->mesh;
    met  = parmesh->met;

    // Make a copy of the mesh before the remeshing
    if ( !PMMG2D_copy_mesh_met ( parmesh ) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to copy old mesh\n",__func__);
      return 0;
    }

    // Reset the value of the fem mode
    mesh->info.fem = parmesh->info.fem;

    // Mark reinitialisation in order to be able to remesh all the mesh
    mesh->mark = 0;
    mesh->base = 0;

    for ( k = 1 ; k <= mesh->nt ; k++ ) {
      mesh->tria[k].flag = mesh->base;
    }

    // Scale the mesh before the remeshing
    if (mesh->nt) {
      ier = MMG5_scaleMesh(mesh,met,NULL);
    }
    else {
      ier = 1;
    }

    MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

    if ( !ieresult ) PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);

    // Remesh with MMG2D
    if (mesh->nt) {
      ier = MMG2D_mmg2d1n(mesh,met);
    }
    else {
      ier = 1;
    }

    mesh->edge = NULL;
    mesh->npi = mesh->np;
    mesh->nti = mesh->nt;

    // Delete the adjacency table and build the new table in MMG2D_pack
    if ( mesh->adja ) PMMG2D_DEL_MEM(mesh,mesh->adja,int,"adja table");

    if (mesh->nt && !MMG2D_pack(mesh,met,NULL)) {
      fprintf(stderr,"\n  ## Triangle packing problem. Exit program.\n");
      ier = 0;
    }

    MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
    if ( !ieresult ) {
      fprintf(stderr,"\n  ## MMG remeshing problem. Exit program.\n");
      PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
    }

    if ( parmesh->iter < parmesh->niter-1 && (!parmesh->info.inputMet) ) {
      // Delete the metric computed by Mmg except at last iter
      PMMG2D_DEL_MEM(mesh,met->m,double,"internal metric");
    }

    // Unscale the mesh
    if (mesh->nt) {
      ier = MMG5_unscaleMesh(mesh,met,NULL) ;
    }
    else {
      ier = 1;
    }

    // MMG2D remove the MG_REQ on boundary nodes
    // Impose again the tags on these nodes
    if (mesh->info.nosurf) {
      for ( k = 1; k <= mesh->nt; k++ ) {

        pt = &mesh->tria[k];

        for ( int i = 0; i <= 2; i++ ) {
          if ( pt->tag[i] & MG_BDY ) {
            pt->tag[i] |= MG_NOSURF;
            pt->tag[i] |= MG_REQ;
          }
        }
      }

      for ( k = 1; k <= mesh->np; k++ ) {
        if ( mesh->point[k].tag & MG_BDY ) {
          mesh->point[k].tag |= MG_NOSURF;
          mesh->point[k].tag |= MG_REQ;
        }
      }
    }

    MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
    if ( !ieresult ) {
      fprintf(stderr,"\n  ## MMG unscaling problem. Exit program.\n");
      PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
    }

    // Copy the metrics on the required edges
    if ( parmesh->iter < parmesh->niter-1 ) {

      if (mesh->nt) {
        np = parmesh->old_mesh->np+5;
      }
      else {
        np = 0;
      }

      MPI_Allreduce(&np, &np_tot, 1, MPI_INT, MPI_SUM, parmesh->comm);

      if (mesh->nt) {
        ier = PMMG2D_copyMetrics_point( mesh, parmesh->old_mesh, met, parmesh->old_met, np_tot );
      }
      else {
        ier = 1;
      }

      MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
      if ( !ieresult ) {
        fprintf(stderr,"\n  ## Copy metrics on required edges problem. Exit program.\n");
        PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
      }
    }

    // Reset the mesh->gap field in case Mmg have modified it
    mesh->gap = MMG5_GAP;

    if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
      chrono(OFF,&(ctim[tim]));
      printim(ctim[tim].gdif,stim);
      fprintf(stdout,"\n       mmg                               %s\n",stim);
    }

    // Interpolate metrics 
    if ( parmesh->iter < parmesh->niter-1 ) {
      if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
        tim = 2;
        chrono(RESET,&(ctim[tim]));
        chrono(ON,&(ctim[tim]));
      }

      if (mesh->nt) {
        ier = PMMG2D_interpMetrics( parmesh, polygon, size );
      }
      else {
        ier = 1;
      }

      MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
      if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
        chrono(OFF,&(ctim[tim]));
        printim(ctim[tim].gdif,stim);
        fprintf(stdout,"       metric interpolation   %s\n",stim);
      }

      if ( !ieresult ) {
        if ( !parmesh->myrank )
          fprintf(stderr,"\n  ## Metrics or fields interpolation problem. Try to save the mesh and exit program.\n");
        PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
      }
    }

    // Free old mesh
    MMG2D_Free_all( MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &parmesh->old_mesh,
                    MMG5_ARG_ppMet,  &parmesh->old_met,
                    MMG5_ARG_end );

    // Renumber the global ref number of vertices so they are unique
    PMMG2D_renumber_mesh( parmesh );

    // Update the list of interface nodes
    PMMG2D_free_interface_nodes_list( parmesh );
    PMMG2D_fill_interface_nodes_list( parmesh );

    if ( parmesh->iter == parmesh->niter-1 ) break;

    // Front advancing
    tim = 3;
    if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
      chrono(RESET,&(ctim[tim]));
      chrono(ON,&(ctim[tim]));
    }

    if (parmesh->nprocs > 1) {

      ier = PMMG2D_frontadvancing( parmesh );

      MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
      if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
        chrono(OFF,&(ctim[tim]));
        printim(ctim[tim].gdif,stim);
        fprintf(stdout,"       front advancing                    %s\n",stim);
      }

      if ( !ieresult ) {
        if ( !parmesh->myrank )
          fprintf(stderr,"\n  ## Front advancing problem. Exit program.\n");
        PMMG2D_CLEAN_AND_RETURN(parmesh, PMMG2D_STRONGFAILURE);
      }

    }

  }

  if (polygon != NULL) {
    for (k = 0; k < size; k++) free(polygon[k]);
    free(polygon);
  }

  // Write a vtk file containing the mesh and the metric
  //check_vtk_file(parmesh);

  if ( parmesh->info.imprim > PMMG2D_VERB_STEPS ) {
    printf("\n");
  }

  // Compute statistical properties of the mesh
  ier = PMMG2D_qualhisto(parmesh, 0);

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    ier_end = PMMG2D_LOWFAILURE;
  }

  if ( parmesh->info.imprim > PMMG2D_VERB_STEPS ) {
    tim = 4;
    chrono(ON,&(ctim[tim]));
  }

  PMMG2D_CLEAN_AND_RETURN(parmesh, ier_end);

}

