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
 * \file libparmmg2d.c
 * \brief Wrapper for the parallel remeshing library.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Main library functions (parallel remeshing starting from centralized data).
 *
 */

#include "parmmg2d.h"
#include "git_log_pmmg2d.h"

// Declared in the header, but defined at compile time 
extern int (*PMMG2D_interp3bar)(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol oldMet, MMG5_pTria ptr, int, PMMG2D_barycoord *barycoord);
extern int (*PMMG2D_interp2bar)(MMG5_pMesh mesh, MMG5_pSol met, MMG5_pSol oldMet, MMG5_pTria ptr, int ip, int l, PMMG2D_barycoord *barycoord);

/**
 * \param parmesh pointer toward the parmesh
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if merge or boundary
 * reconstruction fail
 *
 * Mesh post-treatment:
 *  1. merge meshes if centralized output is asked;
 *  2. build boundary entities of merged or distributed meshes (triangles, edges...);
 *
 */
static inline
int PMMG2D_parmmglib_post(PMMG2D_pParMesh parmesh) {
  mytime        ctim[TIMEMAX];
  int           ier,iresult;
  int8_t        tim;
  char          stim[32];

  tminit(ctim,TIMEMAX);

  iresult = 1;

  switch ( parmesh->info.fmtout ) {
  case ( PMMG2D_UNSET ):
    // No output
    break;
  case ( MMG5_FMT_VtkPvtu ): case ( PMMG2D_FMT_Distributed ): case ( PMMG2D_FMT_Centralized_Distributed ):
  case ( PMMG2D_FMT_DistributedMeditASCII ): case ( PMMG2D_FMT_DistributedMeditBinary ):
  case ( PMMG2D_FMT_Distributed_VtkVtk ):
    int max_nt, min_nt;
    MPI_Allreduce( &parmesh->mesh->nt, &max_nt, 1, MPI_INT, MPI_MAX, parmesh->comm );
    MPI_Allreduce( &parmesh->mesh->nt, &min_nt, 1, MPI_INT, MPI_MIN, parmesh->comm );
    if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
      fprintf(stdout, "\n       Balancing charges - minimal number of local triangles: %d "
                      "- maximal number of local triangles: %d \n", min_nt, max_nt);
    }

    // Check whether repartitioning is necessary
    if (max_nt > (1+parmesh->info.ratio_load_balance) * min_nt) {

      if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
        chrono(RESET,&(ctim[tim]));
        chrono(ON,&(ctim[tim]));
      }

      if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) fprintf( stdout,"\n   -- PHASE 3 : PARTITION THE FINAL MESH \n");

      int sum_nt;
      MPI_Allreduce( &parmesh->mesh->nt, &sum_nt, 1, MPI_INT, MPI_SUM, parmesh->comm );

      // The decomposition is done using parallel tools
      if (sum_nt > PARALLEL_DECOMPOSITION_LIMIT) {

#ifdef USE_PARMETIS
        ier = PMMG2D_ParMetis_decomposition(parmesh);
#elif USE_PTSCOTCH
        ier = PMMG2D_PTScotch_decomposition(parmesh);
#endif

        MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
        if ( !iresult ) {
          if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
#ifdef USE_PARMETIS
            fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO DECOMPOSE THE MESH WITH PARMETIS \n");
#elif USE_PTSCOTCH
            fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO DECOMPOSE THE MESH WITH PT-SCOTCH \n");
#endif
          }
        }

        // Update the list of interface nodes
        PMMG2D_free_interface_nodes_list( parmesh );
        PMMG2D_fill_interface_nodes_list( parmesh );

        MPI_Allreduce( &parmesh->mesh->nt, &max_nt, 1, MPI_INT, MPI_MAX, parmesh->comm );
        MPI_Allreduce( &parmesh->mesh->nt, &min_nt, 1, MPI_INT, MPI_MIN, parmesh->comm );

        if (max_nt > (1+parmesh->info.ratio_load_balance) * min_nt && parmesh->info.imprim > PMMG2D_VERB_VERSION) {
#ifdef USE_PARMETIS
          fprintf(stdout, "\n The ParMetis decomposition does not lead to very balanced partitions - "
                          "minimal number of local triangles: %d - maximal number of local triangles: %d \n", min_nt, max_nt);
#elif USE_PTSCOTCH
          fprintf(stdout, "\n The PT-Scotch decomposition does not lead to very balanced partitions - "
                          "minimal number of local triangles: %d - maximal number of local triangles: %d \n", min_nt, max_nt);
#endif
        }
      }
      // The decomposition is done using sequential tools
      else {

        // Centralize the mesh and then use Scotch or Metis

        ier = PMMG2D_merge_parmesh( parmesh );
        MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

        if ( !iresult ) {
          // Try to save at parallel format
          if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
            fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO CENTRALIZE MESHES..."
                    " TRY TO SAVE DISTRIBUTED MESHES\n\n\n");
          }
        }

        if (parmesh->mesh->adja) PMMG2D_DEL_MEM(parmesh->mesh, parmesh->mesh->adja, int, "remove adjacency table" );

        if (parmesh->myrank) {

          int sol_type = parmesh->met->type;

          MMG2D_Free_all( MMG5_ARG_start,
                          MMG5_ARG_ppMesh, &parmesh->mesh,
                          MMG5_ARG_ppMet,  &parmesh->met,
                          MMG5_ARG_ppDisp, &parmesh->disp,
                          MMG5_ARG_ppLs,   &parmesh->ls,
                          MMG5_ARG_end );

          // Delete the list of interface nodes
          PMMG2D_free_interface_nodes_list(parmesh);

          parmesh->mesh  = NULL;
          parmesh->met   = NULL;
          parmesh->disp  = NULL;
          parmesh->ls    = NULL;

          MMG2D_Init_mesh(MMG5_ARG_start,
                          MMG5_ARG_ppMesh,&parmesh->mesh,
                          MMG5_ARG_ppMet,&parmesh->met,
                          MMG5_ARG_end);

          MMG2D_Set_solSize(parmesh->mesh, parmesh->met, MMG5_Vertex, 0, sol_type);

        }

#ifdef USE_METIS // Use Metis by default as it is far quicker than scotch
        ier = PMMG2D_Metis_decomposition(parmesh);
#elif USE_SCOTCH
        ier = PMMG2D_Scotch_decomposition(parmesh);
#endif

        MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
        if ( !iresult ) {
          if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
#ifdef USE_METIS
            fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO DECOMPOSE THE MESH WITH METIS \n");
#elif USE_SCOTCH
            fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO DECOMPOSE THE MESH WITH SCOTCH \n");
#endif
          }
        }

        // Update the list of interface nodes
        PMMG2D_free_interface_nodes_list( parmesh );
        PMMG2D_fill_interface_nodes_list( parmesh );

        MPI_Allreduce( &parmesh->mesh->nt, &max_nt, 1, MPI_INT, MPI_MAX, parmesh->comm );
        MPI_Allreduce( &parmesh->mesh->nt, &min_nt, 1, MPI_INT, MPI_MIN, parmesh->comm );

        if (max_nt > (1+parmesh->info.ratio_load_balance) * min_nt && parmesh->info.imprim > PMMG2D_VERB_VERSION) {
#ifdef USE_METIS
          fprintf(stdout, "\n The Metis decomposition does not lead to very balanced partitions - "
                          "minimal number of local triangles: %d - maximal number of local triangles: %d \n", min_nt, max_nt);
#elif USE_SCOTCH
          fprintf(stdout, "\n The Scotch decomposition does not lead to very balanced partitions - "
                          "minimal number of local triangles: %d - maximal number of local triangles: %d \n", min_nt, max_nt);
#endif
        }
      }

      if ( parmesh->info.imprim > PMMG2D_VERB_ITWAVES ) {
        chrono(OFF,&(ctim[tim]));
        printim(ctim[tim].gdif,stim);
        fprintf(stdout,"\n       Parallel domain decomposition                  %s\n",stim);
      }

      // If a centralized mesh is also needed
      if (parmesh->info.fmtout == PMMG2D_FMT_Centralized_Distributed) {

        ier = PMMG2D_merge_parmesh_initial( parmesh );

        MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
        if ( !iresult ) {
          if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
            fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO MERGE THE MESH \n");
          }
        }
      }

    }
    // No need of repartitioning but centralized output demanded
    else if (parmesh->info.fmtout == PMMG2D_FMT_Centralized_Distributed) {
      ier = PMMG2D_merge_parmesh_initial( parmesh );
      MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
      if ( !iresult ) {
        if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
          fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO MERGE THE MESH \n");
        }
      }
    }

    break;
  default:
    // Centralized Output
    // Merge all the meshes on the proc 0
    tim = 1;
    chrono(ON,&(ctim[tim]));
    if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
      fprintf( stdout,"\n   -- PHASE 3 : MERGE MESHES OVER PROCESSORS\n" );
    }

    if (parmesh->nprocs > 1) ier = PMMG2D_merge_parmesh( parmesh );
    MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

    if ( !iresult && parmesh->nprocs > 1) {
      // Try to save at parallel format
      if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
        fprintf(stdout,"\n\n\n  -- IMPOSSIBLE TO CENTRALIZE MESHES..."
                " TRY TO SAVE DISTRIBUTED MESHES\n\n\n");
      }
    }

    chrono(OFF,&(ctim[tim]));
    if ( parmesh->info.imprim >  PMMG2D_VERB_VERSION  ) {
      printim(ctim[tim].gdif,stim);
      fprintf( stdout,"   -- PHASE 3 COMPLETED.     %s\n",stim );
    }

  }

  return PMMG2D_SUCCESS;
}

/**
 * \param parmesh pointer toward the parmesh
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if not
 *
 * Distribute the mesh over the processors
 *
 */
int PMMG2D_distributeMesh_centralized_timers( PMMG2D_pParMesh parmesh, mytime *ctim ) {
  MMG5_pMesh    mesh;
  MMG5_pSol     met;
  int           ier,iresult;
  int8_t        tim;
  char          stim[32];

  /** Check input data */
  tim = 1;
  chrono(ON,&(ctim[tim]));

  ier = PMMG2D_check_inputData( parmesh );
  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

  if ( !iresult ) return PMMG2D_LOWFAILURE;

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED.     %s\n",stim);
  }

  chrono(ON,&(ctim[2]));
  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS AND MESH DISTRIBUTION\n");
  }

  /** Mesh preprocessing: set function pointers, scale mesh, perform mesh
   * analysis and display length and quality histos. */
  if( parmesh->myrank == parmesh->info.root ) {
    tim = 7;
    if ( parmesh->info.imprim >= PMMG2D_VERB_STEPS ) {
      chrono(ON,&(ctim[tim]));
      fprintf(stdout,"\n  -- ANALYSIS" );
    }

    ier = PMMG2D_preprocessMesh( parmesh, 1 );

    if ( parmesh->info.imprim >= PMMG2D_VERB_STEPS ) {
      chrono(OFF,&(ctim[tim]));
      printim(ctim[tim].gdif,stim);
      fprintf(stdout,"\n  -- ANALYSIS COMPLETED    %s\n",stim );
    }

    mesh = parmesh->mesh;
    met  = parmesh->met;

    if ( (ier==PMMG2D_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met, NULL ) ) ier = PMMG2D_LOWFAILURE;

    /* Memory repartition */
    if ( !PMMG2D_updateMeshSize( parmesh,1 ) ) ier = 3;

    for (int k = 1; k <= parmesh->mesh->nt; k++) parmesh->mesh->tria[k].ref = k;

  } 
  else {
    ier = PMMG2D_SUCCESS;
  }

  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );

  if ( iresult != PMMG2D_SUCCESS ) {
    return iresult;
  }

  // Send mesh partionning to other procs
  tim = 8;
  if ( parmesh->info.imprim >= PMMG2D_VERB_STEPS ) {
    chrono(ON,&(ctim[tim]));
    fprintf(stdout,"\n  -- PARTITIONING" );
  }

  if ( !PMMG2D_distribute_mesh( parmesh ) ) {
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_LOWFAILURE);
  }

  if ( parmesh->info.imprim >= PMMG2D_VERB_STEPS ) {
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"\n  -- PARTITIONING COMPLETED    %s\n",stim );
  }

  // Function setters (must be assigned before quality computation)
  if( parmesh->myrank != parmesh->info.root ) {
    MMG2D_Set_commonFunc();
    MMG2D_setfunc(parmesh->mesh, parmesh->met);
    PMMG2D_setfunc(parmesh);
  }

  chrono(OFF,&(ctim[2]));
  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    printim(ctim[2].gdif,stim);
    fprintf(stdout,"  -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  iresult = PMMG2D_SUCCESS;
  return iresult;

}

/**
 * \param parmesh pointer toward the parmesh
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if not
 *
 * Perform the mesh remeshing from a centralized mesh
 *
 */
int PMMG2D_parmmg2dlib_centralized(PMMG2D_pParMesh parmesh, double* velocity) {
  int           ier;
  int           ierlib = 0;
  mytime        ctim[TIMEMAX];
  int8_t        tim;
  char          stim[32];

  int remeshing = parmesh->niter;

  if ( parmesh->info.imprim > PMMG2D_VERB_NO ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMG2DLIB_CENTRALIZED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG2D_STR,PMMG2D_VERSION_RELEASE,PMMG2D_RELEASE_DATE,PMMG2D_STR);
    fprintf(stdout,"     git branch: %s\n",PMMG2D_GIT_BRANCH);
    fprintf(stdout,"     git commit: %s\n",PMMG2D_GIT_COMMIT);
    fprintf(stdout,"     git date:   %s\n\n",PMMG2D_GIT_DATE);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  if (!remeshing && parmesh->myrank == parmesh->info.root) {
    for (int k = 1; k <= parmesh->mesh->np; k++) parmesh->mesh->point[k].ref = k;
  }

  // Distribute the mesh
  ier = PMMG2D_distributeMesh_centralized_timers( parmesh, ctim );

  if( ier != PMMG2D_SUCCESS ) return ier;

  // Remeshing
  if (remeshing) {
    tim = 3;
    chrono(ON,&(ctim[tim]));
    if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
      fprintf( stdout,"\n  -- PHASE 2 : %s MESHING\n",
               parmesh->met->size < 3 ? "ISOTROPIC" : "ANISOTROPIC" );
    }
  
    ier = PMMG2D_parmmg2dlib1(parmesh, velocity);
    MPI_Allreduce( &ier, &ierlib, 1, MPI_INT, MPI_MAX, parmesh->comm );
  
    chrono(OFF,&(ctim[tim]));
    printim(ctim[tim].gdif,stim);
    if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
      fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
    }
    if ( ierlib == PMMG2D_STRONGFAILURE ) {
      return ierlib;
    }
  }
  else {
    // Update the list of interface nodes and return
    PMMG2D_free_interface_nodes_list( parmesh );
    PMMG2D_fill_interface_nodes_list( parmesh );
    PMMG2D_CLEAN_AND_RETURN(parmesh, ier);
  }

  ier = PMMG2D_parmmglib_post(parmesh);
  ierlib = MG_MAX ( ier, ierlib );

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"\n   PARMMG2DLIB_CENTRALIZED: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMG2DLIB_CENTRALIZED \n  %s\n",
            PMMG2D_STR, PMMG2D_STR);
  }

  PMMG2D_CLEAN_AND_RETURN(parmesh,ierlib);

}

/**
 * \param parmesh pointer toward the parmesh
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if not
 *
 * Perform the mesh remeshing from a distributed mesh
 *
 */
int PMMG2D_parmmg2dlib_distributed(PMMG2D_pParMesh parmesh, double* velocity) {
  MMG5_pMesh       mesh;
  MMG5_pSol        met;
  int              ier,iresult,ierlib;
  mytime           ctim[TIMEMAX];
  int8_t           tim;
  char             stim[32];

  if ( parmesh->info.imprim >= PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"\n  %s\n   MODULE PARMMG2DLIB_DISTRIBUTED: IMB-LJLL : "
            "%s (%s)\n  %s\n",PMMG2D_STR,PMMG2D_VERSION_RELEASE,PMMG2D_RELEASE_DATE,PMMG2D_STR);
    fprintf(stdout,"     git branch: %s\n",PMMG2D_GIT_BRANCH);
    fprintf(stdout,"     git commit: %s\n",PMMG2D_GIT_COMMIT);
    fprintf(stdout,"     git date:   %s\n\n",PMMG2D_GIT_DATE);
  }

  tminit(ctim,TIMEMAX);
  chrono(ON,&(ctim[0]));

  // Check input data 
  tim = 1;
  chrono(ON,&(ctim[tim]));

  if ( parmesh->info.fmtout == PMMG2D_FMT_Unknown ) {
    parmesh->info.fmtout = PMMG2D_FMT_Centralized;
  }

  ier = PMMG2D_check_inputData( parmesh );
  MPI_CHECK( MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MIN, parmesh->comm ),
             return PMMG2D_LOWFAILURE);
  if ( !iresult ) return PMMG2D_LOWFAILURE;

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"  -- CHECK INPUT DATA COMPLETED.     %s\n",stim);
  }

  tim = 2;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"\n  -- PHASE 1 : ANALYSIS\n");
  }

  // Mesh preprocessing: set function pointers, scale mesh, perform mesh
  // analysis and display length and quality histos.
  ier  = PMMG2D_preprocessMesh( parmesh, 0 );

  mesh = parmesh->mesh;
  met  = parmesh->met;
  if ( (ier == PMMG2D_STRONGFAILURE) && MMG5_unscaleMesh( mesh, met, NULL ) ) {
    ier = PMMG2D_LOWFAILURE;
  }

  MPI_Allreduce( &ier, &iresult, 1, MPI_INT, MPI_MAX, parmesh->comm );
  if ( iresult != PMMG2D_SUCCESS ) {
    return iresult;
  }

  // Track the parallel interface and associate them between processors
  if (parmesh->nprocs > 1) ier = PMMG2D_bind_parallel_interfaces( parmesh );

  chrono(OFF,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    printim(ctim[tim].gdif,stim);
    fprintf(stdout,"   -- PHASE 1 COMPLETED.     %s\n",stim);
  }

  // Remeshing
  tim = 3;
  chrono(ON,&(ctim[tim]));
  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    fprintf( stdout,"\n  -- PHASE 2 : %s MESHING\n",
             met->size < 3 ? "ISOTROPIC" : "ANISOTROPIC" );
  }

  ier = PMMG2D_parmmg2dlib1(parmesh, velocity);
  MPI_Allreduce( &ier, &ierlib, 1, MPI_INT, MPI_MAX, parmesh->comm );

  chrono(OFF,&(ctim[tim]));
  printim(ctim[tim].gdif,stim);
  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"  -- PHASE 2 COMPLETED.     %s\n",stim);
  }
  if ( ierlib == PMMG2D_STRONGFAILURE ) {
    return ierlib;
  }

  // Post-processing (mostly merge the mesh on rank 0)
  ier = PMMG2D_parmmglib_post(parmesh);
  ierlib = MG_MAX ( ier, ierlib );

  chrono(OFF,&ctim[0]);
  printim(ctim[0].gdif,stim);
  if ( parmesh->info.imprim >= PMMG2D_VERB_VERSION ) {
    fprintf(stdout,"\n   PARMMG2DLIB_DISTRIBUTED: ELAPSED TIME  %s\n",stim);
    fprintf(stdout,"\n  %s\n   END OF MODULE PARMMG2DLIB_DISTRIBUTED: IMB-LJLL \n  %s\n",
            PMMG2D_STR, PMMG2D_STR);
  }

  PMMG2D_CLEAN_AND_RETURN(parmesh,ierlib);

}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Check the validity of the input mesh data (tetra orientation, solution
 * compatibility with respect to the provided mesh, Mmg options).
 *
 */
int PMMG2D_check_inputData(PMMG2D_pParMesh parmesh)
{
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  if ( parmesh->info.imprim > PMMG2D_VERB_VERSION )
    fprintf(stdout,"\n  -- PMMG: CHECK INPUT DATA\n");

  mesh = parmesh->mesh;
  met  = parmesh->met;

  // Check options
  if ( mesh->info.lag > -1 ) {
    fprintf(stderr, "  ## Error: lagrangian mode unavailable (MMG2D_IPARAM_lag):\n");
    return 0;
  } 
  else if ( mesh->info.iso ) {
    fprintf(stderr,"  ## Error: level-set discretisation unavailable (MMG2D_IPARAM_iso):\n");
    return 0;
  } 

  // Specific meshing
  if ( met->np ) {
    if ( mesh->info.optim ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: OPTIM OPTION CAN NOT BE USED"
             " WITH AN INPUT METRIC.\n");
      return 0;
    }

    if ( mesh->info.hsiz>0. ) {
      printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ OPTION CAN NOT BE USED"
             " WITH AN INPUT METRIC.\n");
      return 0;
    }
  }

  if ( mesh->info.optim &&  mesh->info.hsiz > 0. ) {
    printf("\n  ## ERROR: MISMATCH OPTIONS: HSIZ AND OPTIM OPTIONS CAN NOT BE USED"
           " TOGETHER.\n");
    return 0;
  }

  // Load data
  if ( met->np && (met->np != mesh->np) ) {
    fprintf(stdout,"  ## WARNING: WRONG METRIC NUMBER. IGNORED\n");
    MMG5_DEL_MEM(mesh,met->m);
    met->np = 0;
  } 
  else if ( met->size!=1 && met->size!=3 ) {
    fprintf(stderr,"  ## ERROR: WRONG DATA TYPE.\n");
    return 0;
  }

  return 1;
}

/**
 * \param  parmesh pointer to parmesh structure
 *
 * \return PMMG_SUCCESS if success, PMMG_LOWFAILURE if fail and return an
 * unscaled mesh, PMMG_STRONGFAILURE if fail and return a scaled mesh.
 *
 * Mesh preprocessing: set function pointers, scale mesh, perform mesh
 * analysis and display length and quality histos.
 */
int PMMG2D_preprocessMesh( PMMG2D_pParMesh parmesh, int is_Central )
{
  int k,i;
  MMG5_pTria pt;
  MMG5_pMesh mesh;
  MMG5_pSol  met;

  mesh = parmesh->mesh;
  met  = parmesh->met;

  assert ( ( mesh != NULL ) && ( met != NULL ) && "Preprocessing empty args");

  // Function setters (must be assigned before quality computation) 
  MMG2D_Set_commonFunc();

  // Mesh scaling and quality histogram
  if ( !MMG5_scaleMesh(mesh,met,NULL) ) {
    return PMMG2D_LOWFAILURE;
  }
  // Don't reset the hmin value computed when unscaling the mesh
  if ( !parmesh->info.sethmin ) {
    mesh->info.sethmin = 1;
  }
  // Don't reset the hmax value computed when unscaling the mesh
  if ( !parmesh->info.sethmax ) {
    mesh->info.sethmax = 1;
  }

  // specific meshing 
  if ( mesh->info.optim && !met->np ) {
    if ( !MMG2D_doSol(mesh,met) ) {
      return PMMG2D_STRONGFAILURE;
    }
  }

  if ( mesh->info.hsiz > 0. ) {
    if ( !MMG2D_Set_constantSize(mesh,met) ) {
      return PMMG2D_STRONGFAILURE;
    }
  }

  MMG2D_setfunc(mesh,met);
  PMMG2D_setfunc(parmesh);

  if ( !PMMG2D_qualhisto(parmesh, is_Central) ) {
    return PMMG2D_STRONGFAILURE;
  }

  for (k=1; k<=mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) )  continue;

    for (i=0; i<3; i++) pt->edg[i] = -1;
  }

  // Mesh analysis
  if ( !MMG2D_analys(mesh) ) {
    return PMMG2D_STRONGFAILURE;
  }

  // Correct the NOSURF tags removed by MMG2D_analysis improperly
  if (!is_Central && !PMMG2D_correct_BDY_tags(parmesh)) {
    return PMMG2D_STRONGFAILURE;
  }

  if ( parmesh->info.imprim0 > PMMG2D_VERB_ITWAVES && (!mesh->info.iso) && met->m ) {
    PMMG2D_prilen(parmesh);
  }

  // Mesh unscaling
  if ( !MMG5_unscaleMesh(mesh,met,NULL) ) {
    return PMMG2D_STRONGFAILURE;
  }

  return PMMG2D_SUCCESS;
}

