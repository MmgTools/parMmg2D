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
 * \file parmmg2d.c
 * \brief main file for the parmmg2d application
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "parmmg2d.h"

mytime         PMMG2D_ctim[TIMEMAX];


// Print elapsed time at end of process.
static void PMMG2D_endcod(void) {
  char   stim[32];

  chrono(OFF,&PMMG2D_ctim[0]);
  printim(PMMG2D_ctim[0].gdif,stim);
  fprintf(stdout,"\n   ELAPSED TIME  %s\n",stim);
}

int main(int argc,char *argv[])
{
  PMMG2D_pParMesh parmesh = NULL;
  int           rank;
  int           ier,ieresult,fmtin,fmtout;
  int8_t        tim, distributedInput;
  char          stim[32],*ptr;

  // Shared memory communicator: processes that are on the same node, sharing
  //    local memory and can potentially communicate without using the network
  MPI_Comm comm_shm = 0;
  int      rank_shm = 0;

  // Initializations: MPI, mesh, and memory
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // Print info with rank 0
  if ( !rank ) {
    fprintf(stdout,"  -- PARMMG2D, Release %s (%s) \n",PMMG2D_VERSION_RELEASE, PMMG2D_RELEASE_DATE);
    fprintf(stdout,"     %s\n",PMMG2D_COPYRIGHT);
    fprintf(stdout,"     %s %s\n\n",__DATE__,__TIME__);

    fprintf(stdout,"  -- MMG2D,    Release %s (%s) \n",MMG_VERSION_RELEASE, MMG_RELEASE_DATE);
    fprintf(stdout,"     %s\n",MMG_COPYRIGHT);
  }

  // Print timer at exit
  if ( !rank ) {
    atexit(PMMG2D_endcod);
  }

  MMG2D_Set_commonFunc();
  tminit(PMMG2D_ctim,TIMEMAX);
  chrono(ON,&PMMG2D_ctim[0]);

  // Allocate the main pmmg2d struct and assign default values 
  if ( PMMG2D_Init_parMesh( PMMG2D_ARG_start,
                            PMMG2D_ARG_ppParMesh, &parmesh,
                            PMMG2D_ARG_dim, 2,
                            PMMG2D_ARG_MPIComm, MPI_COMM_WORLD,
                            PMMG2D_ARG_end) != 1 ) {
    MPI_Abort( MPI_COMM_WORLD, PMMG2D_STRONGFAILURE );
    MPI_Finalize();
    return PMMG2D_FAILURE;
  }

  // reset default values for file names 
  if ( MMG2D_Free_names(MMG5_ARG_start,
                        MMG5_ARG_ppMesh, &parmesh->mesh,
                        MMG5_ARG_ppMet,  &parmesh->met,
                        MMG5_ARG_end) != 1 )
    PMMG2D_RETURN_AND_FREE( parmesh, PMMG2D_STRONGFAILURE );

  // Set default metric size
  if ( !MMG2D_Set_solSize(parmesh->mesh, parmesh->met, MMG5_Vertex, 0, MMG5_Scalar) )
      PMMG2D_RETURN_AND_FREE(parmesh, PMMG2D_STRONGFAILURE);

  if ( PMMG2D_parsar( argc, argv, parmesh ) != 1 ) PMMG2D_RETURN_AND_FREE( parmesh, PMMG2D_STRONGFAILURE );

  if ( parmesh->ddebug ) {
    MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                         &comm_shm );
    MPI_Comm_rank( comm_shm, &rank_shm );
    MPI_Comm_size( comm_shm, &parmesh->size_shm );

    if ( !rank_shm )
      printf("\n     %d MPI PROCESSES (%d ON LOCAL NODE):\n",parmesh->nprocs,
             parmesh->size_shm);

    printf("         MPI RANK %d (LOCAL RANK %d)\n", parmesh->myrank,rank_shm );
  }

  // load data 
  tim = 1;
  chrono(ON,&PMMG2D_ctim[tim]);
  if ( rank==parmesh->info.root && parmesh->info.imprim > PMMG2D_VERB_NO ) {
    fprintf(stdout,"\n  -- INPUT DATA: LOADING MESH ON RANK %d\n",
            parmesh->info.root);
  }

  ptr   = MMG5_Get_filenameExt(parmesh->meshin);

  fmtin = MMG5_Get_format(ptr,MMG5_FMT_MeditASCII);

  // Compute default output format 
  ptr = MMG5_Get_filenameExt(parmesh->meshout);

  // Format from output mesh name 
  fmtout = MMG5_Get_format(ptr,fmtin);

  distributedInput = 0;

  // Load the mesh
  ier = PMMG2D_loadMesh_centralized(parmesh, parmesh->meshin, fmtin);
  MPI_Bcast( &ier, 1, MPI_INT, parmesh->info.root, parmesh->comm );

  if ( 1 != ier ) {
    // try to load distributed mesh
    ieresult = PMMG2D_loadMesh_distributed(parmesh,parmesh->meshin);
    MPI_Allreduce( &ieresult, &ier, 1, MPI_INT, MPI_MIN, parmesh->comm);
    if ( ier == 1) distributedInput = 1;
  }

  if ( ier < 1 ) {
    if ( ier == 0 && rank == parmesh->info.root) {
      fprintf(stderr,"  ** %s  NOT FOUND.\n",parmesh->meshin);
      fprintf(stderr,"  ** UNABLE TO OPEN INPUT FILE.\n");
    }
    PMMG2D_RETURN_AND_FREE(parmesh,PMMG2D_STRONGFAILURE);
  }

  if ( parmesh->mesh->info.lag >= 0 ) {
    if ( rank == parmesh->info.root ) {
      fprintf(stderr,"\n  ## ERROR: LAGRANGIAN MODE I/O NOT YET IMPLEMENTED\n");
    }
    PMMG2D_RETURN_AND_FREE(parmesh, PMMG2D_LOWFAILURE );
  }

  if ( parmesh->mesh->info.lag >= 0 || parmesh->mesh->info.iso ) {
    // displacement or isovalue are mandatory (but not implemented)
    if ( rank == parmesh->info.root ) {
      printf("  ## Error: displacement input not yet implemented.\n");
    }
    PMMG2D_RETURN_AND_FREE(parmesh, PMMG2D_LOWFAILURE );
  }
  else {
    // Facultative metric
    if ( !distributedInput ) {
      ieresult = PMMG2D_loadMet_centralized( parmesh, parmesh->metin );
    }
    else {
      int ier = PMMG2D_loadMet_distributed( parmesh, parmesh->metin );
      MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm);
    }
    if ( ieresult == -1 ) {
      if ( rank == parmesh->info.root ) {
        fprintf(stderr,"\n  ## ERROR: UNABLE TO LOAD METRIC.\n");
      }
      PMMG2D_RETURN_AND_FREE(parmesh, PMMG2D_LOWFAILURE );
    }
  }

  // In iso mode: read metric if any
  if ( parmesh->mesh->info.iso && parmesh->metin ) {
    ier = PMMG2D_loadMet_centralized( parmesh, parmesh->metin );
    if ( ier == -1 ) {
      if ( rank == parmesh->info.root ) {
        fprintf(stderr,"\n  ## ERROR: UNABLE TO LOAD METRIC.\n");
      }
      PMMG2D_RETURN_AND_FREE(parmesh, PMMG2D_LOWFAILURE );
    }
  }

  if ( !PMMG2D_parsop(parmesh) ) PMMG2D_RETURN_AND_FREE(parmesh, PMMG2D_LOWFAILURE );

  chrono(OFF,&PMMG2D_ctim[tim]);
  if ( parmesh->info.imprim > PMMG2D_VERB_NO ) {
    printim(PMMG2D_ctim[tim].gdif,stim);
    fprintf(stdout,"  -- DATA READING COMPLETED.     %s\n",stim);
  }

  // Main call to the remeshing code
  if ( parmesh->mesh->mark ) {
    /* Save a local parameters file containing the default parameters */
    printf("  ## Error: default parameter file saving not yet implemented.\n");
    ier = 2;
    PMMG2D_RETURN_AND_FREE(parmesh,ier);
  }
  else if ( distributedInput ) {
    // Parallel remeshing starting from a distributed mesh
    ier = PMMG2D_parmmg2dlib_distributed(parmesh);
  }
  else {
    // Parallel remeshing starting from a centralized mesh
    ier = PMMG2D_parmmg2dlib_centralized(parmesh);
  }

  /** Check result and save output files */
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( ieresult == PMMG2D_STRONGFAILURE ) { PMMG2D_RETURN_AND_FREE( parmesh, ier ); }

  tim = 2;
  chrono(ON,&PMMG2D_ctim[tim]);

  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"\n  -- WRITING DATA FILE %s\n",parmesh->meshout);
  }

  switch ( parmesh->info.fmtout ) {
  case ( PMMG2D_UNSET ):
    // No output 
    if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
      printf("     ... SKIPPING!\n");
    }
    break;
  case ( MMG5_FMT_VtkPvtu ):
  case ( MMG5_FMT_GmshASCII ): 
  case ( MMG5_FMT_GmshBinary ):
  case ( MMG5_FMT_VtkVtu ):
  case ( MMG5_FMT_VtkVtk ):
    parmesh->meshout[strlen(parmesh->meshout)-6] = '\0';
    strcat(parmesh->meshout, "vtk");
    ier = PMMG2D_saveMeshandMet_centralized_vtk(parmesh);
    if ( !ier ) PMMG2D_RETURN_AND_FREE(parmesh,PMMG2D_STRONGFAILURE);
    break;
  case ( PMMG2D_FMT_Distributed ):
    parmesh->meshout[strlen(parmesh->meshout)-7] = '\0';
    char suffix_file[20];
    sprintf(suffix_file, "%d.mesh", parmesh->myrank);
    strcat(parmesh->meshout, suffix_file);
    parmesh->metout[strlen(parmesh->metout)-6] = '\0';
    sprintf(suffix_file, "%d.sol", parmesh->myrank);
    strcat(parmesh->metout, suffix_file);
    ier = PMMG2D_saveMeshandMet_distributed(parmesh);
    if ( !ier ) PMMG2D_RETURN_AND_FREE(parmesh,PMMG2D_STRONGFAILURE);
    break;
  case ( PMMG2D_FMT_Distributed_VtkVtk ):
    parmesh->meshout[strlen(parmesh->meshout)-7] = '\0';
    sprintf(suffix_file, "%d.vtk", parmesh->myrank);
    strcat(parmesh->meshout, suffix_file);
    ier = PMMG2D_saveMeshandMet_distributed_vtk(parmesh);
    if ( !ier ) PMMG2D_RETURN_AND_FREE(parmesh,PMMG2D_STRONGFAILURE);
    break;
  default:
    ier = PMMG2D_saveMesh_centralized(parmesh,parmesh->meshout);
    if ( !ier ) PMMG2D_RETURN_AND_FREE(parmesh,PMMG2D_STRONGFAILURE);

    if ( parmesh->met && parmesh->met->m ) {
      if ( !PMMG2D_saveMet_centralized(parmesh,parmesh->metout) ) {
        PMMG2D_RETURN_AND_FREE(parmesh,PMMG2D_STRONGFAILURE);
      }
    }

    break;
  }

  MPI_Finalize();

  return 0;
} 
