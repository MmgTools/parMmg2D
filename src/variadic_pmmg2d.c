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
 * \file variadic_pmmg2d.c
 * \brief C variadic functions definitions for PMMG2D library.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \date 03 2024
 * \copyright GNU Lesser General Public License.
 *
 * \note This file contains some internal functions for the API, see
 * the \ref libparmmg2d.h header file for the documentation of all
 * the usefull user's API functions.
 *
 * variadic functions definitions for PMMG2D library.
 *
 */

#include "parmmg2d.h"

/**
 * \param argptr list of the type of structures that must be initialized inside
 * your parmesh and needed informations for ParMmg2D (mesh dimension and MPI
 * communicator).
 *
 * \param callFromC 1 if called from C API, 0 if called from the Fortran one.
 *
 * \a argptr contains at least a pointer toward a parmesh pointer preceeded by
 * the PMMG2D_ARG_pParMesh keyword and the mesh dimension preceeded by the
 * PMMG2D_ARG_dim keyword
 *
 * By default, ParMmg2D will initilize at least 1 mesh, 1 metric and the
 * MPI_COMM_WORLD_COMMUNICATOR inside your parmesh. Thus, the 2 following calls
 * are identicals:
 *
 * 1) MMG2D_Init_parmesh(PMMG2D_ARG_start,PMMG2D_ARG_ppParMesh,your_pParmesh_address,
 *    PMMG2D_ARG_pMesh,PMMG2D_ARG_pMet,
 *    PMMG2D_ARG_dim,mesh_dimension,PMMG2D_ARG_MPIComm,MPI_COMM_WORLD,PMMG2D_ARG_end)
 *
 * 2) MMG2D_Init_parmesh(PMMG2D_ARG_start,PMMG2D_ARG_ppParMesh,your_pParmesh_address,
 *    PMMG2D_ARG_dim,mesh_dimension,PMMG2D_ARG_end)
 *
 * \return 1 if success, 0 if fail
 *
 * Internal function for structure allocations (taking a va_list argument).
 *
 */

int PMMG2D_Init_parMesh_var_internal(va_list argptr, int callFromC ) {
  PMMG2D_pParMesh  *parmesh;
  MPI_Comm         comm;
  int    typArg,dim,comm_f;
  int    parmeshCount,meshCount,metCount,dimCount,commCount;
  int    lsCount,dispCount;

  parmeshCount = 0;
  meshCount    = 0;
  metCount     = 0;
  lsCount      = 0;
  dispCount    = 0;
  dimCount     = 0;
  commCount    = 0;

  dim  = 2;
  comm = MPI_COMM_WORLD;
  while ( (typArg = va_arg(argptr,int)) != PMMG2D_ARG_end )
  {
    switch ( typArg )
    {
    case(PMMG2D_ARG_ppParMesh):
      parmesh = va_arg(argptr,PMMG2D_pParMesh*);
      ++parmeshCount;
      break;
    case(PMMG2D_ARG_pMesh):
      ++meshCount;
      break;
    case(PMMG2D_ARG_pMet):
      ++metCount;
      break;
    case(PMMG2D_ARG_pLs):
      ++lsCount;
      break;
    case(PMMG2D_ARG_pDisp):
      ++dispCount;
      break;
    case(PMMG2D_ARG_dim):
      ++dimCount;
      dim = va_arg(argptr,int);
      break;
    case(PMMG2D_ARG_MPIComm):
      ++commCount;
      if ( callFromC ) {
        comm = va_arg(argptr,MPI_Comm);
      }
      else {
        comm_f = va_arg(argptr,int);
        comm = MPI_Comm_f2c(comm_f);
      }
      break;
    default:
      fprintf(stderr,"\n  ## Error: PMMG2D_Init_parmesh:\n"
              " unexpected argument type: %s\n",PMMG2D_Get_pmmg2dArgName(typArg));
      return 0;
    }
  }

  if ( parmeshCount != 1 ) {
    fprintf(stderr,"\n  ## Error: PMMG2D_Init_parmesh:\n"
            " you need to initialize the parmesh structure that"
            " will contain your data (mesh, metric, communicator...\n");
    return 0;
  }

  if ( meshCount > 1 ) {
    fprintf(stdout,"\n  ## Warning: PMMG2D_Init_parmesh:\n"
            " Only 1 mesh structure is allowed.\n");
  }
  if ( metCount > 1 ) {
    fprintf(stdout,"\n  ## Warning: PMMG2D_Init_parmesh:\n"
            " Only 1 metric structure is allowed.\n");
  }
  if ( lsCount != 0 ) {
    fprintf(stderr,"\n  ## Error: PMMG2D_Init_parmesh:\n"
            " Level-set is not yet implemented in ParMMG2D.\n");
  }
  if ( dispCount != 0 ) {
    fprintf(stderr,"\n  ## Error: PMMG2D_Init_parmesh:\n"
            " Displacement is not yet implemented in ParMMG2D.\n");
  }
  if ( commCount > 1 ) {
    fprintf(stdout,"\n  ## Warning: PMMG2D_Init_parmesh:\n"
            " More than 1 MPI communicator provided. Used the last one.\n");
  }
  if ( dimCount > 1 ) {
    fprintf(stdout,"\n  ## Warning: PMMG2D_Init_parmesh:\n"
            " More than 1 dimension provided. Used the last one.\n");
  }
  else if ( !dimCount ) {
    fprintf(stderr,"\n  ## Error: PMMG2D_Init_parmesh:\n"
            " you need to provide the dimension of your mesh using the PMMG2D_dim"
            " keyword\n.");
    return 0;
  }

  if ( dim != 2 ) {
    fprintf(stderr,"\n  ## Error: PMMG2D_Init_parmesh:\n"
            " dimension must be equal to 2 with ParMMG2D.\n"
            " use ParMMG3D if the mesh is 3D.\n");
    return 0;
  }

  // ParMesh allocation 
  assert ( (*parmesh == NULL) && "trying to initialize non empty parmesh" );
  MMG5_SAFE_CALLOC( *parmesh, 1, PMMG2D_ParMesh, return 0 );

  if ( *parmesh == NULL ) {
    return 0;
  }

  // Assign some values to memory related fields to begin working with
  (*parmesh)->memGloMax = 4 * 1024L * 1024L;
  (*parmesh)->memMax = 4 * 1024L * 1024L;
  (*parmesh)->memCur = sizeof(PMMG2D_ParMesh);

  (*parmesh)->mesh  = NULL;
  (*parmesh)->met   = NULL;
  (*parmesh)->disp  = NULL;
  (*parmesh)->ls    = NULL;

  // Init the mesh and the metric (no level-set and displacement)
  if ( !MMG2D_Init_mesh(MMG5_ARG_start,
                        MMG5_ARG_ppMesh,&(*parmesh)->mesh,
                        MMG5_ARG_ppMet,&(*parmesh)->met,
                        MMG5_ARG_end) )
      return MMG5_STRONGFAILURE;

  PMMG2D_Init_parameters(*parmesh,comm);

  return 1;
}

/**
 * \param argptr list of the parmmg2d structures that must be deallocated.
 *
 * \a argptr contains at least a pointer toward a \a PMMG2D_pParMesh structure
 * (that will contain the parmesh and identified by the PMMG2D_ARG_ppParMesh keyword)
 *
 * \return 0 if fail, 1 if success
 *
 * Deallocations of the parmmg2d structures before return
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 */
int PMMG2D_Free_all_var(va_list argptr)
{

  PMMG2D_pParMesh  *parmesh;
  int              typArg;
  int              parmeshCount;

  parmeshCount = 0;

  while ( (typArg = va_arg(argptr,int)) != PMMG2D_ARG_end )
  {
    switch ( typArg )
    {
    case(PMMG2D_ARG_ppParMesh):
      parmesh = va_arg(argptr,PMMG2D_pParMesh*);
      ++parmeshCount;
      break;
    default:
      fprintf(stdout,"\n  ## Warning: PMMG2D_Free_all:\n"
              " ignored argument: %s\n",PMMG2D_Get_pmmg2dArgName(typArg));
    }
  }

  PMMG2D_Free_names( *parmesh );

  MMG2D_Free_all( MMG5_ARG_start,
                  MMG5_ARG_ppMesh, &(*parmesh)->mesh,
                  MMG5_ARG_ppMet,  &(*parmesh)->met,
                  MMG5_ARG_ppDisp, &(*parmesh)->disp,
                  MMG5_ARG_ppLs,   &(*parmesh)->ls,   
                  MMG5_ARG_end );

  if ((*parmesh)->old_mesh != NULL) 
    MMG2D_Free_all( MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &(*parmesh)->old_mesh,
                    MMG5_ARG_ppMet,  &(*parmesh)->old_met,
                    MMG5_ARG_end );

  if ((*parmesh)->info.optim_interp && (*parmesh)->initial_mesh != NULL) {
    MMG2D_Free_all( MMG5_ARG_start,
                    MMG5_ARG_ppMesh, &(*parmesh)->initial_mesh,
                    MMG5_ARG_ppMet,  &(*parmesh)->initial_met,
                    MMG5_ARG_end );
  }

  // Delete the list of interface nodes
  PMMG2D_free_interface_nodes_list(*parmesh);

  (*parmesh)->memCur -= sizeof(PMMG2D_ParMesh);

  if ( (*parmesh)->info.imprim>5 || (*parmesh)->ddebug ) {
    printf("  MEMORY USED AT END (Bytes) %zu\n",(*parmesh)->memCur);
  }

  MMG5_SAFE_FREE(*parmesh);

  return 1;

}
