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
 * \file parmmg2d.h
 * \brief internal functions headers for parmmg2d
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef _PARMMG2D_H
#define _PARMMG2D_H

#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <signal.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <mpi_pmmg2d.h>

#include "libparmmg2d.h"
#include "interpmesh_pmmg2d.h"
#include "libmmg2d.h"
#include "mmgcommon_private.h"
#include "libmmg2d_private.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \def PMMG2D_NUL
 *
 * Null value
 *
 */
#define PMMG2D_NUL     0

/**
 * \def PMMG2D_EPS
 *
 * Epsilon value
 *
 */
#define PMMG2D_EPS     1e-10

/**
 * \def PMMG2D_NITER
 *
 * Default number of iterations
 *
 */
#define PMMG2D_NITER  4

/**
 * \def PMMG2D_IMPRIM
 *
 * Default verbosity
 *
 */
#define PMMG2D_IMPRIM   5

/** 
 *
 * Number of elements layers for metric interpolation
 * Equal 0 for interpolation on only the triangle containing the points 
 *
 */
#define PMMG2D_INTERP_NLAYERS  0

/**
 * \def PMMG2D_MMG_IMPRIM
 *
 * Default verbosity for Mmg
 *
 */
#define PMMG2D_MMG_IMPRIM   -1

/**
 * \def PMMG2D_MMG_NOSIZREQ
 *
 * Default value of info.nosizreq for Mmg
 * It is strongly advised to keep the default value
 * of PMMG2D
 *
 */
#define PMMG2D_MMG_NOSIZREQ 1

/**
 * \def PMMG2D_MMG_HGRADREQ
 *
 * Default value of info.hgradreq for Mmg
 * It is strongly advised to keep the default value
 * of PMMG2D
 *
 */
#define PMMG2D_MMG_HGRADREQ -1

/** 
 *
 * Number of elements layers for interface displacement 
 *
 */
#define PMMG2D_MVIFCS_NLAYERS 5

/**
 * \def PMMG2D_OPTIM_INTERP
 *
 * Parameter to use the initial mesh for interpolation
 * May fail with domains containing holes
 *
 */
#define PMMG2D_OPTIM_INTERP 0

/**
 *
 * Size of quality histogram arrays
 *
 */
#define PMMG2D_QUAL_HISSIZE 5

/**
 *
 * Admissible load unbalacing at the end of the remeshing process (in % of variation) 
 *
 */
#define PMMG2D_RATIO_LOAD_BALANCE 0.3

/** 
 *
 * Tolerance to load balance when using Scotch, PT-Scotch or ParMetis for partitioning (in %) 
 *
 */
#define PMMG2D_LOAD_IMBALANCE 0.05

/**
 *
 * no verbosity for pmmg library
 *
 */
#define PMMG2D_VERB_NO -1

/**
 *
 * minimal verbosity for pmmg library: print library version and duration
 *
 */
#define PMMG2D_VERB_VERSION 0

/**
 *
 * average verbosity for pmmg library: add parmmg steps information
 *
 */
#define PMMG2D_VERB_STEPS 3

/**
 *
 * average verbosity for pmmg library: add waves information
 *
 */
#define PMMG2D_VERB_ITWAVES 4

/**
 *
 * detailed verbosity for pmmg library: add detailed quality histo
 *
 */
#define PMMG2D_VERB_DETQUAL 5

/**
 *
 * Use custom partitioning saved in the reference field (1=yes, 0=no)
 *
 */
#define PMMG2D_PREDEF_PART 0

/** 
 *
 * Strategy to partition with PT-Scotch (0: default, 1: HPC, 2: good balance, 3: "quick" calculation)
 *
 */
#define SCOTCH_strategy 3

/**
 *
 * Threshold for the use PT-Scotch or ParMetis instead of Scotch or Metis (number of triangles)
 *
 */
#define PARALLEL_DECOMPOSITION_LIMIT 5e7

/**
 * \enum PMMG2D_Format
 * \brief Type of supported file format
 */
enum PMMG2D_Format {
  PMMG2D_FMT_MeditASCII  = MMG5_FMT_MeditASCII, /*!< ASCII Medit (.mesh) */
  PMMG2D_FMT_MeditBinary = MMG5_FMT_MeditBinary,/*!< Binary Medit (.meshb) */
  PMMG2D_FMT_GmshASCII   = MMG5_FMT_GmshASCII,  /*!< ASCII Gmsh */
  PMMG2D_FMT_GmshBinary  = MMG5_FMT_GmshBinary, /*!< Binary Gmsh */
  PMMG2D_FMT_VtkPvtp     = MMG5_FMT_VtkPvtp,    /*!< VTK pvtp */
  PMMG2D_FMT_VtkPvtu     = MMG5_FMT_VtkPvtu,    /*!< VTK pvtu */
  PMMG2D_FMT_VtkVtu      = MMG5_FMT_VtkVtu,     /*!< VTK vtu */
  PMMG2D_FMT_VtkVtp      = MMG5_FMT_VtkVtp,     /*!< VTK vtp */
  PMMG2D_FMT_VtkVtk      = MMG5_FMT_VtkVtk,     /*!< VTK vtk */
  PMMG2D_FMT_Tetgen      = MMG5_FMT_Tetgen,     /*!< Tetgen or Triangle */
  PMMG2D_FMT_Centralized,                       /*!< Centralized Setters/Getters */
  PMMG2D_FMT_Distributed,                       /*!< Distributed Setters/Getters */
  PMMG2D_FMT_DistributedMeditASCII,             /*!< Distributed ASCII Medit (.mesh) */
  PMMG2D_FMT_DistributedMeditBinary,            /*!< Distributed Binary Medit (.meshb) */
  PMMG2D_FMT_Distributed_VtkVtk,                /*!< Distributed VTK vtk */
  PMMG2D_FMT_Centralized_Distributed,           /*!< Centralized and Distributed Setters/Getters */
  PMMG2D_FMT_Unknown,                           /*!< Unrecognized */
};

/**
 * \param parmesh pointer toward a parmesh structure
 * \param val     exit value
 *
 * Controlled parmmg termination:
 *   Deallocate parmesh struct and its allocated members
 *   If this is an unsuccessful exit call abort to cancel any remaining processes
 *   Call MPI_Finalize / exit
 */

#define PMMG2D_RETURN_AND_FREE(parmesh,val) do                            \
  {                                                                     \
                                                                        \
    if ( !PMMG2D_Free_all( PMMG2D_ARG_start,                                \
                         PMMG2D_ARG_ppParMesh,&parmesh,                   \
                         PMMG2D_ARG_end) ) {                              \
      fprintf(stderr,"  ## Warning: unable to clean the parmmg memory.\n" \
              " Possible memory leak.\n");                              \
    }                                                                   \
                                                                        \
    MPI_Finalize();                                                     \
    return(val);                                                        \
                                                                        \
  } while(0)

/**
 * Clean the mesh, the metric and the solutions and return \a val.
 */
#define PMMG2D_CLEAN_AND_RETURN(parmesh,val)do                            \
  {                                                                     \
      if ( parmesh->mesh ) {                                            \
        parmesh->mesh->npi = parmesh->mesh->np; \
        parmesh->mesh->nti = parmesh->mesh->nt; \
        parmesh->mesh->nai = parmesh->mesh->na; \
        parmesh->mesh->nei = parmesh->mesh->ne; \
      }                                                                 \
                                                                        \
      if ( parmesh->met ) parmesh->met->npi  = parmesh->met->np;        \
                                                                        \
    return val;                                                         \
                                                                        \
  }while(0)


#define ERROR_AT(msg1,msg2)                                          \
  fprintf( stderr, msg1 msg2 " function: %s, file: %s, line: %d \n", \
           __func__, __FILE__, __LINE__ )

#define MEM_CHK_AVAIL(mesh,bytes,msg) do {                            \
  if ( (mesh)->memCur + (bytes) > (mesh)->memMax ) {                  \
    ERROR_AT(msg," Exceeded max memory allowed: ");      \
    stat = PMMG2D_FAILURE;                                              \
  } else if ( (mesh)->memCur + (bytes) < 0  ) {                       \
    ERROR_AT(msg," Tried to free more mem than allocated: " );        \
    stat = PMMG2D_SUCCESS;                                              \
  }                                                                   \
  else {                                                              \
    stat = PMMG2D_SUCCESS;                                              \
  } } while(0)

#define PMMG2D_DEL_MEM(mesh,ptr,type,msg) do {                \
    size_t size_to_free;                                    \
                                                            \
    if ( ptr ) {                                            \
      size_to_free = myfree( ptr );                         \
      assert ( (mesh)->memCur >= size_to_free );            \
      (mesh)->memCur -= size_to_free;                       \
      (ptr) = NULL;                                         \
    }                                                       \
  } while(0)

#define PMMG2D_MALLOC(mesh,ptr,size,type,msg,on_failure) do { \
  int    stat = PMMG2D_SUCCESS;                               \
  size_t size_to_allocate;                                  \
                                                            \
  (ptr) = NULL;                                             \
  if ( (size) != 0 ) {                                      \
    size_to_allocate = (size)*sizeof(type);                 \
    MEM_CHK_AVAIL(mesh,size_to_allocate,msg );              \
    if ( stat == PMMG2D_SUCCESS ) {                           \
      (ptr) = (type*)mymalloc( size_to_allocate );          \
      if ( (ptr) == NULL ) {                                \
        ERROR_AT( msg, " malloc failed: " );                \
        on_failure;                                         \
      } else {                                              \
        (mesh)->memCur += size_to_allocate;                 \
        stat = PMMG2D_SUCCESS;                                \
      }                                                     \
    } else {                                                \
      on_failure;                                           \
    }                                                       \
  } } while(0)

#define PMMG2D_CALLOC(mesh,ptr,size,type,msg,on_failure) do { \
  int    stat = PMMG2D_SUCCESS;                               \
  size_t size_to_allocate;                                  \
                                                            \
  (ptr) = NULL;                                             \
  if ( (size) != 0 ) {                                      \
    size_to_allocate = (size)*sizeof(type);                 \
    MEM_CHK_AVAIL(mesh,size_to_allocate,msg);               \
    if ( stat == PMMG2D_SUCCESS ) {                           \
      (ptr) = (type*)mycalloc( (size), sizeof(type) );      \
      if ( (ptr) == NULL ) {                                \
        ERROR_AT(msg," calloc failed: ");                   \
        on_failure;                                         \
      } else {                                              \
        (mesh)->memCur += size_to_allocate;                 \
      }                                                     \
    } else {                                                \
      on_failure;                                           \
    }                                                       \
  } } while(0)

#define PMMG2D_REALLOC(mesh,ptr,newsize,oldsize,type,msg,on_failure) do { \
  int    stat = PMMG2D_SUCCESS;                                           \
  size_t size_to_allocate,size_to_add,size_to_increase;                 \
  type*  tmp;                                                           \
                                                                        \
  if ( (ptr) == NULL ) {                                                \
    assert(((oldsize)==0) && "NULL pointer pointing to non 0 sized memory?"); \
    PMMG2D_MALLOC(mesh,ptr,(newsize),type,msg,on_failure);                \
  } else if ((newsize)==0) {                                            \
    PMMG2D_DEL_MEM(mesh,ptr,type,msg);                                    \
  } else if ((newsize) < (oldsize)) {                                   \
    size_to_allocate = (newsize)*sizeof(type);                          \
    tmp = (type *)myrealloc((ptr),size_to_allocate,                     \
                            (oldsize)*sizeof(type));                    \
    if ( tmp == NULL ) {                                                \
      ERROR_AT(msg," Realloc failed: ");                                \
      PMMG2D_DEL_MEM(mesh,ptr,type,msg);                                  \
      on_failure;                                                       \
    } else {                                                            \
      (ptr) = tmp;                                                      \
      (mesh)->memCur -= (((oldsize)*sizeof(type))-size_to_allocate);    \
    }                                                                   \
  } else if ((newsize) > (oldsize)) {                                   \
    size_to_add = ((newsize)-(oldsize))*sizeof(type);                   \
    size_to_allocate = (newsize)*sizeof(type);                          \
    size_to_increase = (oldsize)*sizeof(type);                          \
                                                                        \
    MEM_CHK_AVAIL(mesh,size_to_add,msg);                                \
    if ( stat == PMMG2D_SUCCESS ) {                                       \
      tmp = (type *)myrealloc((ptr),size_to_allocate,size_to_increase); \
      if ( tmp == NULL ) {                                              \
        ERROR_AT(msg, " Realloc failed: " );                            \
        PMMG2D_DEL_MEM(mesh,ptr,type,msg);                                \
        on_failure;                                                     \
      } else {                                                          \
        (ptr) = tmp;                                                    \
        (mesh)->memCur += ( size_to_add );                              \
      }                                                                 \
    }                                                                   \
    else {                                                              \
      on_failure;                                                       \
    }                                                                   \
  }                                                                     \
  } while(0)

#define PMMG2D_RECALLOC(mesh,ptr,newsize,oldsize,type,msg,on_failure) do { \
    int my_stat = PMMG2D_SUCCESS;                                         \
                                                                        \
    PMMG2D_REALLOC(mesh,ptr,newsize,oldsize,type,msg,my_stat=PMMG2D_FAILURE;on_failure;); \
    if ( (my_stat == PMMG2D_SUCCESS ) && ((newsize) > (oldsize)) ) {      \
      memset( (ptr) + oldsize, 0, ((size_t)((newsize)-(oldsize)))*sizeof(type)); \
    }                                                                   \
  } while(0)


/**
 * \param parmesh pointer toward a parmesh structure
 * \param mesh pointer toward a mesh structure
 * \param on_failure instruction to execute if fail
 *
 * Check the allowed memory. */
#define PMMG2D_MEM_CHECK(parmesh,mesh,on_failure) do {             \
    size_t memGloMax,memMax;                                      \
    memGloMax = parmesh->memGloMax;                               \
    memMax = mesh->memMax;                                        \
    if( memMax != memGloMax ) {                                   \
      fprintf(stderr,"\n  ## Error: %s: allowed memory mismatch." \
                     " Maximal: %zu -- global: %zu\n",            \
              __func__,memMax,memGloMax);                         \
      on_failure;                                                 \
    }                                                             \
  } while(0)


// Input
int PMMG2D_Set_name(PMMG2D_pParMesh,char **,const char* name,const char* defname);
int PMMG2D_check_inputData ( PMMG2D_pParMesh parmesh );
int PMMG2D_preprocessMesh( PMMG2D_pParMesh parmesh, int isCentral );
int PMMG2D_parsar( int argc, char *argv[], PMMG2D_pParMesh parmesh );

// Internal library
void PMMG2D_setfunc( PMMG2D_pParMesh parmesh );
int PMMG2D_parmmg2dlib1 ( PMMG2D_pParMesh parmesh, double* velocity );

// Parallel Interface
void PMMG2D_fill_interface_nodes_list( PMMG2D_pParMesh parmesh );
void PMMG2D_free_interface_nodes_list( PMMG2D_pParMesh parmesh );
int PMMG2D_mark_parallel_interface_nodes( PMMG2D_pParMesh parmesh );
int PMMG2D_correct_BDY_tags( PMMG2D_pParMesh parmesh );

// Mesh interpolation
int PMMG2D_interpMetrics( PMMG2D_pParMesh parmesh, double** polygon, int size );
int PMMG2D_copyMetrics_point( MMG5_pMesh mesh, MMG5_pMesh oldMesh, MMG5_pSol met, MMG5_pSol oldMet, int npt);
double** build_partition_contour(MMG5_pMesh mesh, int *size);
int*** grid_size_triangles(MMG5_pMesh mesh, double minX, double minY, double maxX, double maxY, int GRID_SIZE);
int PMMG2D_interpFields( PMMG2D_pParMesh parmesh, double* field );

// Move interfaces
int PMMG2D_frontadvancing( PMMG2D_pParMesh parmesh );
int PMMG2D_remove_points( PMMG2D_pParMesh parmesh );
int PMMG2D_add_triangles( PMMG2D_pParMesh parmesh, MMG5_pTria pt_global,
                          MMG5_pPoint pp_global, double* mm_global, int npt, int npp );
int PMMG2D_check_contiguity( PMMG2D_pParMesh parmesh, int nt_tmp );
void PMMG2D_update( PMMG2D_pParMesh parmesh );
void PMMG2D_untag_parallel( PMMG2D_pParMesh parmesh );
int PMMG2D_exchange( PMMG2D_pParMesh parmesh, int size, int* list_index_pt, MMG5_Tria** list_pt,
                     int* list_index_pp, MMG5_Point** list_pp, double** list_mm,
                     MMG5_pTria* pt_global, MMG5_pPoint* pp_global, double** mm_global, int* npt_global, int* npp_global );

// Memory
int  PMMG2D_parmesh_SetMemMax( PMMG2D_pParMesh parmesh);
int  PMMG2D_setMeshSize( MMG5_pMesh,int,int,int,int );
int  PMMG2D_setMeshSize_alloc( MMG5_pMesh );
int  PMMG2D_setMeshSize_realloc( MMG5_pMesh,int,int,int,int);
int  PMMG2D_updateMeshSize( PMMG2D_pParMesh parmesh,int fitMesh);
void PMMG2D_parmesh_SetMemGloMax( PMMG2D_pParMesh parmesh );

// Tools 
int PMMG2D_copy_mmgInfo ( MMG5_Info *info, MMG5_Info *info_cpy );
int PMMG2D_copy_initial_mesh_met( PMMG2D_pParMesh parmesh );

// Quality
int PMMG2D_qualhisto( PMMG2D_pParMesh parmesh, int );
int PMMG2D_prilen( PMMG2D_pParMesh parmesh );

// Variadic_pmmg2d.c
int PMMG2D_Init_parMesh_var_internal(va_list argptr,int callFromC);
int PMMG2D_Free_all_var(va_list argptr);

const char* PMMG2D_Get_pmmg2dArgName(int typArg);

// Mesh decomposition
#ifdef USE_SCOTCH
int PMMG2D_Scotch_decomposition( PMMG2D_pParMesh parmesh );
int PMMG2D_Scotch_decomposition_root( PMMG2D_pParMesh parmesh );
#endif
#ifdef USE_PTSCOTCH
int PMMG2D_PTScotch_decomposition( PMMG2D_pParMesh parmesh );
int PMMG2D_PTScotch_exchange( PMMG2D_pParMesh parmesh, SCOTCH_Num *partTab );
#endif
int PMMG2D_exchange_from_root( PMMG2D_pParMesh parmesh, int size, MMG5_pTria* pt_global, MMG5_pPoint* pp_global, 
            			       double** mm_global, int* npt_global, int* npp_global );
int PMMG2D_Metis_exchange(PMMG2D_pParMesh parmesh);
int PMMG2D_Metis_decomposition(PMMG2D_pParMesh parmesh);
int PMMG2D_ParMetis_decomposition( PMMG2D_pParMesh parmesh );
int PMMG2D_exchange_mesh( PMMG2D_pParMesh parmesh );
int PMMG2D_find_neighbour_triangles( PMMG2D_pParMesh parmesh, int* nb_neighbours, int **list_neighbours );

#ifdef __cplusplus
}
#endif

#endif
