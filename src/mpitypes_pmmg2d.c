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
 * \file mpitypes_pmmg2d.c
 * \brief Mpi types management.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 */

#include "parmmg2d.h"
#include "mpitypes_pmmg2d.h"

/**
 * \param mpi_point new MPI data type
 *
 * Create an MPI data type named mpi_point to allow the communication of
 * choosen fields of a point (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG2D_create_MPI_Point(MPI_Datatype *mpi_point)
{
  MMG5_Point   point[2];
  int          i,blck_lengths[5] = {3, 3, 1, 1, 1};
  MPI_Aint     displs[5],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT16_T};

  MPI_Get_address(&(point[0]),       &lb);
  MPI_Get_address(&(point[0].c[0]),  &displs[0]);
  MPI_Get_address(&(point[0].n[0]),  &displs[1]);
  MPI_Get_address(&(point[0].ref),   &displs[2]);
  MPI_Get_address(&(point[0].xp),    &displs[3]);
  MPI_Get_address(&(point[0].tag),   &displs[4]);
  MPI_Get_address(&(point[1]),       &ub);

  // Adress relative to the object adress
  for ( i=4 ; i>= 0; --i ) displs[i] -= lb;

  MPI_Type_create_struct(5, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_point);

  MPI_Type_commit(mpi_point);
}

/**
 * \param mpi_edge new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a edge (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG2D_create_MPI_Edge(MPI_Datatype *mpi_edge)
{
  MMG5_Edge    edge[2];
  int          i,blck_lengths[4] = {1, 1, 1, 1};
  MPI_Aint     displs[4],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[4] = {MPI_INT,MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(edge[0]),      &lb);
  MPI_Get_address(&(edge[0].a),    &displs[0]);
  MPI_Get_address(&(edge[0].b),    &displs[1]);
  MPI_Get_address(&(edge[0].ref),  &displs[2]);
  MPI_Get_address(&(edge[0].tag),  &displs[3]);
  MPI_Get_address(&(edge[1]),      &ub);

  // Adress relative to the object adress
  for ( i=3 ; i>= 0; --i ) displs[i] -= lb;

  MPI_Type_create_struct(4, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_edge);

  MPI_Type_commit(mpi_edge);
}

/**
 * \param mpi_tria new MPI data type
 *
 * Create an MPI data type to allow the communication of
 * choosen fields of a tria (only those who needs to be communicated).
 *
 * \warning to fill when we need to communicate additional things
 */
void PMMG2D_create_MPI_Tria(MPI_Datatype *mpi_tria)
{
  MMG5_Tria    tria[2];
  int          i,blck_lengths[3] = {3, 1, 3};
  MPI_Aint     displs[3],lb,ub;
  MPI_Datatype mpi_noextent;
  MPI_Datatype types[3] = {MPI_INT,MPI_INT,MPI_INT16_T};

  MPI_Get_address(&(tria[0]),       &lb);
  MPI_Get_address(&(tria[0].v[0]),  &displs[0]);
  MPI_Get_address(&(tria[0].ref),   &displs[1]);
  MPI_Get_address(&(tria[0].tag[0]),&displs[2]);
  MPI_Get_address(&(tria[1]),       &ub);

  // Adress relative to the object adress
  for ( i=2 ; i>= 0; --i ) displs[i] -= lb;

  MPI_Type_create_struct(3, blck_lengths, displs, types, &mpi_noextent);

  MPI_Type_create_resized(mpi_noextent,lb,ub-lb,mpi_tria);

  MPI_Type_commit(mpi_tria);
}


/**
 * \param mpi_point     pointer toward an MPI_Datatype
 * \param mpi_edge      pointer toward an MPI_Datatype
 * \param mpi_triangle  pointer toward an MPI_Datatype
 *
 * If used, free the \a mpi_point, \a mpi_xpoint, ... MPI_Datatype
 *
 */
void PMMG2D_Free_MPI_meshDatatype( MPI_Datatype *mpi_point,
                                   MPI_Datatype *mpi_edge,
                                   MPI_Datatype *mpi_triangle,
                                   bool edge ) {

  if ( *mpi_triangle ) {
    MPI_Type_free( mpi_triangle );
    *mpi_triangle = MPI_DATATYPE_NULL;
  }

  if ( edge ) {
    MPI_Type_free( mpi_edge );
    *mpi_edge = MPI_DATATYPE_NULL;
  }

  if ( *mpi_point ) {
    MPI_Type_free( mpi_point );
    *mpi_point = MPI_DATATYPE_NULL;
  }

}
