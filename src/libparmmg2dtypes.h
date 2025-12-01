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
 * \file libparmmgtypes.h
 * \brief parmmg types and functions that must be accessible to the library users
 * \author Cécile Dobrzynski (Bx INP/Inria)
 * \author Algiane Froehly (Inria)
 * \version 5
 * \copyright GNU Lesser General Public License.
 */
#ifndef _LIBPARMMG2DTYPES_H
#define _LIBPARMMG2DTYPES_H

#include "mmg/common/libmmgtypes.h"
#include "pmmg2dversion.h"
#include <mpi.h>


/**
 * \def PMMG2D_SUCCESS
 *
 * Return value for success.
 */
#define PMMG2D_SUCCESS       0
/**
 * \def PMMG2D_LOWFAILURE
 *
 * Return value if the remesh process failed but we can save a conform
 * mesh.
 */
#define PMMG2D_LOWFAILURE    1
/**
 * \def PMMG2D_STRONGFAILURE
 *
 * Return value if the remesh process failed and the mesh is
 * non-conform.
 */
#define PMMG2D_STRONGFAILURE 2
/**
 * \def PMMG2D_FAILURE
 *
 * Return value of failure, caller is to decide how to proceed further,
 * ie whether remesh process failed and/or the mesh is non-conformant
 */
#define PMMG2D_FAILURE  4


/**
 * \def PMMG2D_ARG_start
 *
 * To begin a list of variadic arguments (mandatory first arg for all our
 * variadic functions)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_start  1
/**
 * \def PMMG2D_ARG_ppParMesh
 *
 * pointer toward a pointer toward a parMesh structure
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_ppParMesh  2
/**
 * \def PMMG2D_ARG_pMesh
 *
 * PMMG2D_pMesh structure
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_pMesh  3
/**
 * \def PMMG2D_ARG_pMet
 *
 * PMMG2D_pSol structure storing a metric field
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_pMet   4
/**
 * \def PMMG2D_ARG_pSols
 *
 * PMMG2D_pSol structure storing an array of solutions
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_pSols  5
/**
 * \def PMMG2D_ARG_pDisp
 *
 * PMMG2D_pSol structure storing a displacement field
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_pDisp  6
/**
 * \def PMMG2D_ARG_pLs
 *
 * PMMG2D_pSol structure storing level-set function
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_pLs  7
/**
 * \def PMMG2D_ARG_ngroups
 *
 * Number of groups per processor.
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_MPIComm  8
/**
 * \def PMMG2D_ARG_dim
 *
 * mesh dimension
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_dim  9
/**
 * \def PMMG2D_ARG_end
 *
 * To end a list of variadic argument (mandatory last argument for all our
 * variadic functions)
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define PMMG2D_ARG_end    10
/**
 * \def PMMG2D_LOADBALANCING_metis
 *
 * Use metis to compute and balance the graph during the loadbalancing step
 *
 */
#define PMMG2D_LOADBALANCING_metis 1
/**
 * \def PMMG2D_UNSET
 *
 * Initialization value
 *
 */
#define PMMG2D_UNSET     -1

/**
 * Types
 */
/**
 * \struct PMMG2D_int_nodes
 * \brief interface nodes between processors
 */
typedef struct {
  int* index; /*!< index of the node in the other partitions */ 
  int* proc; /*!< processor number associated to the index */
  double* coord; /*!< coordinates of the vertex */
  int nb; /*!< size */
} PMMG2D_int_nodes;
typedef PMMG2D_int_nodes * PMMG2D_pint_nodes;

/**
 * \struct PMMG2D_Info
 * \brief Store input parameters of the run.
 */
typedef struct {

  int imprim;  /*!< ParMmg verbosity (may be non-null only on zero rank) */
  int imprim0; /*!< ParMmg verbosity of the zero rank */
  int mem;     /*!< memory asked by user */
  int iso;     /*!< ls mode (not yet available) */
  int root;    /*!< MPI root rank */
  int fem;     /*!< fem mesh (no elt with more than 1 bdy face */
  int mmg_imprim; /*!< 1 if the user has manually setted the mmg verbosity */
  int ifc_layers;  /*!< nb of layers for interface displacement */
  int optim_interp; /*!< force to use the initial mesh during metric interpolation when possible */
  int interp_layers; /*!< nb of layers for interpolation */
  double ratio_load_balance; /*!< limit for final load balancing */
  int nobalancing; /*!< switch off final load balancing */
  int loadbalancing_mode; /*!< way to perform the loadbalanding (see LOADBALANCING) */
  int contiguous_mode; /*!< force/don't force partitions contiguity */
  int fmtout; /*!< store the output format asked */
  int8_t sethmin; /*!< 1 if user set hmin, 0 otherwise (needed for multiple library calls) */
  int8_t sethmax; /*!< 1 if user set hmin, 0 otherwise (needed for multiple library calls) */
  uint8_t inputMet; /* 1 if User prescribe a metric or a size law */
} PMMG2D_Info;


/**
 * \struct PMMG2D_ParMesh
 * \brief ParMmg mesh structure.
 */
typedef struct {

  /* mpi info */
  MPI_Comm    comm;   /*!< Global communicator of all parmmg processes */
  int         nprocs; /*!< Number of processes in global communicator */
  int         myrank; /*!< Rank in global communicator */
  int         size_shm; /*!< Number or MPI process per Node */
  int         nip; // Number of parallel interface nodes

  /* mem info */
  size_t    memGloMax; /*!< Maximum memory available to all structs */
  size_t    memMax; /*!< Maximum memory parmesh is allowed to allocate */
  size_t    memCur; /*!< Currently allocated memory */

  /* file names */
  char     *meshin,*meshout;
  char     *metin,*metout;
  char     *lsin;
  char     *dispin;

  MMG5_pMesh   mesh;  /*!< mesh definition : coordinates, triangles, etc.. */
  MMG5_pMesh   old_mesh;
  MMG5_pMesh   initial_mesh;
  MMG5_pSol    met;   /*!< metric */
  MMG5_pSol    old_met;
  MMG5_pSol    initial_met;
  MMG5_pSol    disp;  /*!< displacement */ // TODO: to implement
  MMG5_pSol    ls;    /*!< level-set */ // TODO: to implement

  PMMG2D_pint_nodes   list_interface_nodes;

  /* global variables */
  int            ddebug; //! Debug level
  int            iter;   //! Current adaptation iteration
  int            niter;  //! Number of adaptation iterations

  /* parameters of the run */
  PMMG2D_Info      info; /*!< \ref PMMG2D_Info structure */

} PMMG2D_ParMesh;
typedef PMMG2D_ParMesh  * PMMG2D_pParMesh;

#endif
