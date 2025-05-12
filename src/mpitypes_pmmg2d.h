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

#ifndef MPITYPES_PMMG2D_H
#define MPITYPES_PMMG2D_H
/**
 * \file mpitypes_pmmg2d.h
 * \brief Mpi types management header file.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 */
#include <mpi_pmmg2d.h>
#include "libmmgtypes.h"

void PMMG2D_create_MPI_Point(MPI_Datatype *mpi_point);

void PMMG2D_create_MPI_Edge(MPI_Datatype *mpi_edge);

void PMMG2D_create_MPI_Tria(MPI_Datatype *mpi_tria);

void PMMG2D_Free_MPI_meshDatatype( MPI_Datatype*,MPI_Datatype*,MPI_Datatype*,bool);

#endif
