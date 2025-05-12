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
 * \file locate_pmmg2d.h
 * \brief Point localization for interpolation on a new mesh.
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef LOCATE_PMMG2D_H

#define LOCATE_PMMG2D_H

#include "barycoord_pmmg2d.h"

int PMMG2D_precompute_nodeTrias( PMMG2D_pParMesh parmesh,MMG5_pMesh mesh,int **nodeTrias );
int PMMG2D_locatePointInTria( MMG5_pMesh mesh,MMG5_pTria ptr,int k,MMG5_pPoint ppt,
                              PMMG2D_barycoord *barycoord,
                              double *closestDist,int *closestTria );
int PMMG2D_locatePoint( MMG5_pMesh mesh, int *list_trianles, MMG5_pPoint ppt,
                        PMMG2D_barycoord *barycoord,
                        int *iTria,int *foundWedge,int *foundCone );
void PMMG2D_locatePoint_errorCheck( MMG5_pMesh mesh,int ip,int ier,int myrank );
void PMMG2D_locate_setStart( MMG5_pMesh mesh,MMG5_pMesh meshOld );

#endif
