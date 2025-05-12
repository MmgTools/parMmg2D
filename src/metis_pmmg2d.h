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
 * \file metis_pmmg2d.h
 * \brief metis_pmmg2d.c header file
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef METIS_PMMG2D_H

#define METIS_PMMG2D_H

#include "parmmg2d.h"

#ifdef USE_METIS
#include <metis.h>
#endif

#ifdef USE_PARMETIS
#include <parmetis.h>
#endif

/**
 * \def PMMG2D_WGTVAL_HUGEINT
 *
 * Huge integer weight for parallel faces in load balancing
 *
 */
#define PMMG2D_WGTVAL_HUGEINT   1000000

/**
 * \def PMMG2D_CONTIG_DEF
 *
 * value for option[METIS_OPTION_CONTIG] in Metis
 *
 */
#define PMMG2D_CONTIG_DEF     1

int PMMG2D_graph_meshElts2metis(PMMG2D_pParMesh,MMG5_pMesh,MMG5_pSol,idx_t**,idx_t**,idx_t**,idx_t*);
int PMMG2D_part_meshElts2metis(PMMG2D_pParMesh,idx_t*,idx_t);

#endif
