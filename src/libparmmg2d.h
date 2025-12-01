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
 * \file libparmmg.h
 * \brief API headers for the parmmg2d library
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef _PMMG2DLIB_H
#define _PMMG2DLIB_H

#include "libparmmg2dtypes.h"
#ifdef USE_METIS
#include "metis.h"
#endif
#ifdef USE_PARMETIS
#include "parmetis.h"
#endif
#ifdef USE_SCOTCH
#include "scotch.h"
#endif
#ifdef USE_PTSCOTCH
#include "ptscotch.h"
#endif

#if defined(c_plusplus) || defined(__cplusplus)
extern "C" {
#endif

#define PMMG2D_STR   "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"

/**
 * \enum PMMG2D_Param
 * \brief Input parameters for mmg library.
 *
 * Input parameters for mmg library. Options prefixed by \a
 * PMMG2D_IPARAM asked for integers values ans options prefixed by \a
 * PMMG2D_DPARAM asked for real values.
 *
 */
enum PMMG2D_Param {
  PMMG2D_IPARAM_verbose,           /*!< [-10..10], Tune level of verbosity */
  PMMG2D_IPARAM_mmgVerbose,        /*!< [-10..10], Tune level of verbosity of Mmg */
  PMMG2D_IPARAM_mem,               /*!< [n/-1], Set memory size to n Mbytes or keep the default value */
  PMMG2D_IPARAM_debug,             /*!< [1/0], Turn on/off debug mode */
  PMMG2D_IPARAM_distributedOutput, /*!< [0/1], Turn off/on distributed output */
  PMMG2D_IPARAM_mmgDebug,          /*!< [1/0], Turn on/off debug mode */
  PMMG2D_IPARAM_angle,             /*!< [1/0], Turn on/off angle detection */
  PMMG2D_IPARAM_iso,               /*!< [1/0], Level-set meshing */
  PMMG2D_IPARAM_interp_layers,     /*!< [n],  Number of layers for interpolation */
  PMMG2D_IPARAM_lag,               /*!< [-1/0/1/2], Lagrangian option */
  PMMG2D_IPARAM_opnbdy,            /*!< [0/1], Enable preservation of open boundaries */
  PMMG2D_IPARAM_optim,             /*!< [1/0], Optimize mesh keeping its initial edge sizes */
  PMMG2D_IPARAM_nofem,             /*!< [1/0], Generate a non finite element mesh */
  PMMG2D_IPARAM_noinsert,          /*!< [1/0], Avoid/allow point insertion */
  PMMG2D_IPARAM_noswap,            /*!< [1/0], Avoid/allow edge or face flipping */
  PMMG2D_IPARAM_nomove,            /*!< [1/0], Avoid/allow point relocation */
  PMMG2D_IPARAM_nosurf,            /*!< [1/0], Avoid/allow surface modifications */
  PMMG2D_IPARAM_numberOfLocalParam,/*!< [n], Number of local parameters */
  PMMG2D_IPARAM_anisosize,         /*!< [1/0], Turn on/off anisotropic metric creation when no metric is provided */
  PMMG2D_IPARAM_nobalancing,       /*!< [1/0], Deactivate load balancing of the output mesh */
  PMMG2D_IPARAM_ifcLayers,         /*!< [n], Number of layers of interface displacement */
  PMMG2D_IPARAM_niter,             /*!< [n], Set the number of remeshing iterations */
  PMMG2D_DPARAM_angleDetection,    /*!< [val], Value for angle detection */
  PMMG2D_DPARAM_load_balance,      /*!< [val], Value for the tolerance of final load balancing */
  PMMG2D_DPARAM_hmin,              /*!< [val], Minimal mesh size */
  PMMG2D_DPARAM_hmax,              /*!< [val], Maximal mesh size */
  PMMG2D_DPARAM_hsiz,              /*!< [val], Constant mesh size */
  PMMG2D_DPARAM_hausd,             /*!< [val], Control global Hausdorff distance (on all the boundary surfaces of the mesh) */
  PMMG2D_DPARAM_hgrad,             /*!< [val], Control gradation */
  PMMG2D_DPARAM_hgradreq,          /*!< [val], Control gradation from required entities */
  PMMG2D_DPARAM_ls,                /*!< [val], Value of level-set */
};


/* API_functions_pmmg.c */
/* init structures */
/**
 * \param starter dummy argument used to initialize the variadic argument list
 * \param ... variadic arguments that depend to the parmesh fields that you want
 * to init
 *
 * You need to provide at least the following arguments:
 * the \a PMMG2D_ARG_start keyword to start the list of variadic arguments
 * the \a PMMG2D_ARG_ppParMesh keyword to say that the next argument is a pointer
 *  toward a pointer toward a parmesh
 * a pointer toward a pointer toward a parmesh
 * the \a PMMG2D_ARG_pMesh keyword to initialize a \a mesh pointer inside your \a parmesh
 * the \a PMMG2D_ARG_pMet keyword to initialize a \a metric pointer inside your \a parmesh
 * the \a PMMG2D_ARG_dim keyword to set the mesh dimension
 * the \a PMMG2D_ARG_MPIComm keyword to set the MPI Communicator in which parmmg will work
 * the \a PMMG2D_ARG_end keyword to end the list of variadic args.
 *
 * Example :
 * PMMG2D_Init_parmesh(PMMG2D_ARG_start,PMMG2D_ARG_ppParMesh,your_parmesh_address,
 *   PMMG2D_ARG_pMesh,PMMG2D_ARG_pMet,PMMG2D_ARG_dim,mesh_dimension,
 *   PMMG2D_ARG_MPIComm,mpi_communicator,PMMG2D_ARG_end)
 *
 * \return 1 if success, 0 if fail
 *
 * ParMMG structures allocation and initialization.
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 */
 int PMMG2D_Init_parMesh(const int starter,...);

/* libparmmg.c */
/**
 * \param parmesh pointer toward the parmesh structure (boundary entities are
 * stored into MMG5_Tria, MMG5_Edge... structures)
 *
 * \return \ref PMMG2D_SUCCESS if success, \ref PMMG2D_LOWFAILURE if fail but we can
 * one unscaled mesh per proc or \ref PMMG2D_STRONGFAILURE if fail and
 * we can't return one unscaled mesh per proc.
 *
 * Main program for the parallel remesh library for distributed meshes.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_parmmglib_distributed(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 **/
int PMMG2D_parmmg2dlib_distributed(PMMG2D_pParMesh parmesh, double* velocity);

/**
 * \param parmesh pointer toward the parmesh structure (boundary entities are
 * stored into MMG5_Tria, MMG5_Edge... structures)
 *
 * \return \ref PMMG2D_SUCCESS if success, \ref PMMG2D_LOWFAILURE if fail but we can
 * return a centralized and unscaled mesh or \ref PMMG2D_STRONGFAILURE if fail and
 * we can't return a centralized and unscaled mesh.
 *
 * Main program for the parallel remesh library for centralized meshes
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_parmmg2dlib_centralized(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 **/
int PMMG2D_parmmg2dlib_centralized(PMMG2D_pParMesh parmesh, double* velocity);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param comm MPI communicator for ParMmg
 *
 * Initialization of the input parameters (stored in the Info structure).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_INIT_PARAMETERS(parmesh,comm)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER,INTENT(IN)            :: comm\n
 * >   END SUBROUTINE\n
 *
 */
void  PMMG2D_Init_parameters(PMMG2D_pParMesh parmesh,MPI_Comm comm);

/**
 * \param parmesh   Pointer towards the parmesh structure.
 * \param np        Number of vertices.
 * \param ne        Number of tetrahedra.
 * \param nprism    Number of prisms.
 * \param nt        Number of triangles.
 * \param nquad     Number of quads.
 * \param na        Number of edges.
 * \return          0 if failed, 1 otherwise.
 *
 * Set the number of vertices, tetrahedra, prisms, triangles, quadrilaterals and
 * edges of the mesh and allocate the associated tables.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SET_MESHSIZE(parmesh, np, ne, nprism, nt, nquad, na, &\n
 * >                                retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: np, ne, nprism, nt, nquad, na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_Set_meshSize(PMMG2D_pParMesh parmesh, int np, int ne, int nprism, int nt,
                      int nquad, int na);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param iparam integer parameter to set (see \a MMG3D_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SET_IPARAMETER(parmesh,iparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: iparam,val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG2D_Set_iparameter(PMMG2D_pParMesh parmesh, int iparam, int val);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param dparam double parameter to set (see \a MMG3D_Param structure).
 * \param val value for the parameter.
 * \return 0 if failed, 1 otherwise.
 *
 * Set double parameter \a dparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SET_DPARAMETER(parmesh,dparam,val,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(IN)           :: dparam\n
 * >     REAL(KIND=8), INTENT(IN)      :: val\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int  PMMG2D_Set_dparameter(PMMG2D_pParMesh parmesh, int iparam, double val);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Free names stored in the parmesh
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_FREE_NAMES(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_Free_names(PMMG2D_pParMesh parmesh);

/**
 * \param starter dummy argument used to initialize the variadic argument list.
 * \param ... variadic arguments to list the structure that must be deallocated.
 * For now, you must provide the PMMG2D_ARG_ppParMesh keyword and your parmesh
 * address as arguments.
 *
 * \return 1 if success, 0 if fail
 *
 * Deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
  int PMMG2D_Free_all(const int starter,...);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Distribute the mesh from a centralized mesh on root process.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_DISTRIBUTEMESH_CENTRALIZED(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_distributeMesh_centralized(PMMG2D_pParMesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Distribute the mesh from a centralized mesh on root process.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_DISTRIBUTE_MESH(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_distribute_mesh(PMMG2D_pParMesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Attribute unique reference number to each vertex
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_RENUMBER_MESH(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_renumber_mesh(PMMG2D_pParMesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param list_interface_nodes pointer toward the interface nodes
 * \param size_interface number of interface nodes
 * \return 0 if failed, 1 otherwise.
 *
 * Attribute unique reference number to each vertex from distributed meshes
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_RENUMBER_DISTRIBUTED_MESH(parmesh,list_interface_nodes,size,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     MMG5_DATA_PTR_T,INTENT(IN)    :: list_interface_nodes\n
 * >     INTEGER, INTENT(IN)           :: size\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_renumber_distributed_mesh( PMMG2D_pParMesh parmesh, 
                                      PMMG2D_pint_nodes list_interface_nodes, 
                                      int size_interface );

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Track and find the parallel interfaces between the distributed meshes.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_BIND_PARALLEL_INTERFACES(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_bind_parallel_interfaces ( PMMG2D_pParMesh parmesh );

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Set integer parameter \a iparam at value \a val.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_MERGE_PARMESH(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_merge_parmesh( PMMG2D_pParMesh parmesh );
int PMMG2D_merge_parmesh_initial( PMMG2D_pParMesh parmesh );

/**
 * \param parmesh  pointer toward the group structure.
 * \param ghosts table of the neighbouring partitions of each triangle.
 * \return 0   if failed, 1 otherwise.
 *
 * Get ghost elements.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_GET_GHOSTS(parmesh,ghosts,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)           :: parmesh\n
 * >     INTEGER, DIMENSION(*), INTENT(OUT)      :: ghosts\n
 * >     INTEGER, INTENT(OUT)                    :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_Get_ghosts(PMMG2D_pParMesh parmesh, int** ghosts, int* size);

/* libparmmg_tools.c: Tools for the library */
/**
 * \param parmesh pointer to pmmg structure
 * \return 0 if fail, 1 if success.
 *
 * Print the default parameters values
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_DEFAULTVALUES(parmesh,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_defaultValues( PMMG2D_pParMesh parmesh );

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 1 on success
 *         0 on failure
 *
 * Parse command line arguments.
 *
 * \remark no matching fortran function.
 * \remark each proc read the parameters.
 *
 */
int PMMG2D_parsar( int argc, char *argv[], PMMG2D_pParMesh parmesh );

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 1 on success
 *         0 on failure
 *
 * Parse parameter file.
 *
 * \remark no matching fortran function.
 * \remark each proc read the file of parameters.
 *
 */
int PMMG2D_parsop( PMMG2D_pParMesh parmesh );

/**
 * \param parmesh pointer toward the parmesh structure
 * \param prog pointer toward the program name.
 * \param return 1 if success, 0 if fail.
 *
 * Print help for parmmg options.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_USAGE(parmesh,prog,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: prog\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_usage( PMMG2D_pParMesh parmesh, char * const prog);

/* input/output functions */
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 1 if success, 0 or -1 otherwise
 *
 * Read distributed mesh data. Insert rank index to the mesh name.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_LOADMESH_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_loadMesh_distributed(PMMG2D_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 1 if success, 0 or -1 otherwise
 *
 * Read centralized mesh data.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_LOADMESH_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_loadMesh_centralized(PMMG2D_pParMesh parmesh,const char *filename, int ftmin);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load metric field. The solution file must contains only 1 solution: the
 * metric
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_LOADMET_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_loadMet_centralized(PMMG2D_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load metric field. The solution file must contains only 1 solution: the
 * metric. Insert rank index in the file name.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_LOADMET_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_loadMet_distributed(PMMG2D_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load displacement field. The solution file must contains only 1 solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_LOADDISP_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_loadDisp_centralized(PMMG2D_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load level-set field. The solution file must contains only 1 solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_LOADLS_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_loadLs_centralized(PMMG2D_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load displacement, level-set or metric field depending on the
 * option setted. The solution file must contains only 1 solution.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_LOADSOL_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_loadSol_centralized(PMMG2D_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return -1 data invalid, 0 no file, 1 ok.
 *
 * Load 1 or more solutions in a solution file at medit file format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_LOADALLSOLS_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_loadAllSols_centralized(PMMG2D_pParMesh parmesh,const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename pointer toward the name of file.

 * \return 0 if failed, 1 otherwise.
 *
 * Save mesh data for a centralized mesh.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SAVEMESH_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_saveMesh_centralized(PMMG2D_pParMesh parmesh, const char *filename);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SAVEMET_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_saveMet_centralized(PMMG2D_pParMesh parmesh, const char *filename);
/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write isotropic or anisotropic metric of a distributed mesh (insert rank
 * index to filename).
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SAVEMET_DISTRIBUTED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_saveMet_distributed(PMMG2D_pParMesh parmesh, const char *filename);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of file.
 * \return 0 if failed, 1 otherwise.
 *
 * Write 1 or more than 1 solution in a file at medit format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SAVEALLSOLS_CENTRALIZED(parmesh,filename,strlen,retval)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >     CHARACTER(LEN=*), INTENT(IN)   :: filename\n
 * >     INTEGER, INTENT(IN)            :: strlen\n
 * >     INTEGER, INTENT(OUT)           :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_saveAllSols_centralized(PMMG2D_pParMesh parmesh, const char *filename);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and isotropic or anisotropic metric in the vtk format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SAVEALLSOLS_CENTRALIZED(parmesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_saveMeshandMet_centralized_vtk(PMMG2D_pParMesh parmesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and isotropic or anisotropic metric in the vtk format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SAVEMESHANDMET_DISTRIBUTED_VTK(parmesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_saveMeshandMet_distributed_vtk(PMMG2D_pParMesh parmesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if failed, 1 otherwise.
 *
 * Write mesh and isotropic or anisotropic metric in the mesh format.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SAVEMESHANDMET_DISTRIBUTED(parmesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT) :: parmesh\n
 * >   END SUBROUTINE\n
 *
 */
  int PMMG2D_saveMeshandMet_distributed(PMMG2D_pParMesh parmesh);

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * Set function pointers.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_SETFUNC(parmesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: parmesh\n
 * >   END SUBROUTINE\n
 *
 */
void PMMG2D_setfunc( PMMG2D_pParMesh parmesh );

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * Set function pointers.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE PMMG2D_PTScotch_decomposition(parmesh)\n
 * >     MMG5_DATA_PTR_T, INTENT(INOUT)     :: parmesh\n
 * >   END SUBROUTINE\n
 *
 */
int PMMG2D_PTScotch_decomposition( PMMG2D_pParMesh parmesh );

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

#endif
