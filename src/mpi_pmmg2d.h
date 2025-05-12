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

#ifndef MPI_PMMG2D_H
#define MPI_PMMG2D_H
/**
 * \file mpi_pmmg2d.h
 * \brief Mpi tools that are used in different part of parMMG2D
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 */
#include <mpi.h>

#define MPI_ANALYS_TAG                 10000

#define MPI_CHECK(func_call,on_failure) do {                            \
    int mpi_ret_val;                                                    \
                                                                        \
    mpi_ret_val = func_call;                                            \
                                                                        \
    switch ( mpi_ret_val ) {                                            \
    case ( MPI_SUCCESS ):                                               \
      break;                                                            \
    case ( MPI_ERR_COMM ):                                              \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI communicator\n.",   \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_TYPE ):                                              \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI datatype\n.",       \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_COUNT ):                                             \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI count argument\n.", \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_TAG ):                                               \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI tag\n.",            \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_RANK ):                                              \
      fprintf(stderr," ## Error: %s:%d: Invalid MPI rank\n.",           \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    case ( MPI_ERR_OTHER ):                                             \
      fprintf(stderr," ## Error: %s:%d: MPI error.\n.",                 \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    default:                                                            \
      fprintf(stderr," ## Error: %s:%d: Unexpected MPI error\n.",       \
              __func__,__LINE__ );                                      \
      on_failure;                                                       \
      break;                                                            \
    }                                                                   \
  } while(0)
#endif

#define RUN_ON_ROOT_AND_BCAST(func_call,root,myrank,on_failure) do {    \
    int ierloc;                                                         \
                                                                        \
    if ( myrank == root ) {                                             \
      ierloc = func_call;                                               \
    }                                                                   \
    MPI_CHECK( MPI_Bcast(&ierloc,1,MPI_INT,root,parmesh->comm),on_failure); \
    if ( !ierloc ) {                                                    \
      on_failure;                                                       \
    }                                                                   \
                                                                        \
  } while(0)
