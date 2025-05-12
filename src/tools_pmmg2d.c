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
 * \file tools_pmmg2d.c
 * \brief Various tools for parMMG2D
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * various tools for parmmg2d.
 *
 */

#include "parmmg2d.h"

/**
 * \param typArg integer defined by the PMMG2D_ARG_* preprocessor commands
 *
 * \return the name of the \a PMMG2D_ARG_* variable
 *
 * Print the \a PMMG2D_ARG_* name associated to the \a typArg value.
 *
 *
 */
const char* PMMG2D_Get_pmmg2dArgName(int typArg)
{

  switch ( typArg )
  {
  case(PMMG2D_ARG_start):
    return "PMMG2D_ARG_start";
    break;

  case(PMMG2D_ARG_ppParMesh):
    return "PMMG2D_ARG_ppParMesh";
    break;

  case(PMMG2D_ARG_pMesh):
    return "PMMG2D_ARG_pMesh";
    break;

  case(PMMG2D_ARG_pMet):
    return "PMMG2D_ARG_pMet";
    break;

  case(PMMG2D_ARG_pSols):
    return "PMMG2D_ARG_pSols";
    break;

  case(PMMG2D_ARG_pDisp):
    return "PMMG2D_ARG_pDisp";
    break;

  case(PMMG2D_ARG_pLs):
    return "PMMG2D_ARG_pLs";
    break;

  case(PMMG2D_ARG_dim):
    return "PMMG2D_ARG_dim";
    break;

  case(PMMG2D_ARG_MPIComm):
    return "PMMG2D_ARG_MPIComm";
    break;

  case(PMMG2D_ARG_end):
    return "PMMG2D_ARG_end";
    break;

  default:
    return "PMMG2D_ARG_Unknown";
  }
}

/**
 * \param info pointer toward a MMG5 info structure that we want to copy
 * \param info_cpy pointer toward a MMG5 info structure into copy the data
 *
 * \return 1 if success, 0 if not
 *
 * Copy the info data into info_cpy, allocate the \a mat and \a par structures
 * if needed.
 *
 */
int PMMG2D_copy_mmgInfo ( MMG5_Info *info, MMG5_Info *info_cpy ) 
{
  MMG5_pMat mat_tmp;
  MMG5_pPar par_tmp;

  // assert to remove (we may authorize to have mat and par already allocated )
  assert ( (!info_cpy->mat) && (!info_cpy->par) );

  if ( info->nmat && (!info_cpy->mat) ) {
    MMG5_SAFE_CALLOC(mat_tmp,info->nmat,MMG5_Mat,return 0);
  }
  else {
    mat_tmp = info_cpy->mat;
  }
  if ( mat_tmp ) {
    *mat_tmp = *info->mat;
  }

  /* local parameters */
  if ( info->npar && !info_cpy->par ) {
    MMG5_SAFE_CALLOC(par_tmp,info->npar,MMG5_Par,return 0);
  }
  else {
    par_tmp = info_cpy->par;
  }
  if ( par_tmp ) {
    *par_tmp = *info->par;
  }

  *info_cpy = *info;

  info_cpy->mat = mat_tmp;
  info_cpy->par = par_tmp;

  return 1;
}
