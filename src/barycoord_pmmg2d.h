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
 * \file barycoord_pmmg2d.h
 * \brief Barycentric coordinates for point localization in a mesh.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef BARYCOORD_PMMG2D_H

#define BARYCOORD_PMMG2D_H

/** \struct PMMG2D_barycoord
 *
 * \brief Struct containing the index and value of a barycentric coordinate
 *
 */
typedef struct {
  int    idx; /*!< direction */
  double val; /*!< coordinate value */
} PMMG2D_barycoord;


double PMMG2D_quickarea(double *a,double *b,double *c);

void PMMG2D_barycoord_get( double *val, PMMG2D_barycoord *phi, int ndim );

int PMMG2D_barycoord_compare( const void *a,const void *b );

int PMMG2D_barycoord2d_compute( MMG5_pMesh mesh, MMG5_pTria ptr, double *coord,
                                PMMG2D_barycoord *barycoord );

int PMMG2D_barycoord2d_getClosest( MMG5_pMesh mesh,int k,MMG5_pPoint ppt,
                                   PMMG2D_barycoord *barycoord );

int PMMG2D_barycoord_isBorder( PMMG2D_barycoord *phi,int *ifoundEdge,int *ifoundVertex );

int PMMG2D_barycoord2d_evaluate( MMG5_pMesh mesh,MMG5_pTria ptr, double *coord,
                                 PMMG2D_barycoord *barycoord );

#endif
