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
 * \file API_functions_pmmg.c
 * \brief C API functions definitions for PARMMG library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG library.
 *
 */
#include "parmmg2d.h"
#include "metis_pmmg2d.h"
#include "linkedlist_pmmg2d.h"

int PMMG2D_Init_parMesh(const int starter,...) 
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = PMMG2D_Init_parMesh_var_internal(argptr,1);

  va_end(argptr);

  return ier;
}

int PMMG2D_Free_names(PMMG2D_pParMesh parmesh)
{
  PMMG2D_DEL_MEM ( parmesh, parmesh->meshin,char,"meshin" );
  PMMG2D_DEL_MEM ( parmesh, parmesh->meshout,char,"meshout" );
  PMMG2D_DEL_MEM ( parmesh, parmesh->metin,char,"metin" );
  PMMG2D_DEL_MEM ( parmesh, parmesh->metout,char,"metout" );
  PMMG2D_DEL_MEM ( parmesh, parmesh->lsin,char,"lsin" );
  PMMG2D_DEL_MEM ( parmesh, parmesh->dispin,char,"dispin" );
  return 1;
}

int PMMG2D_Free_all(const int starter,...)
{
  va_list argptr;
  int     ier;

  va_start(argptr, starter);

  ier = PMMG2D_Free_all_var(argptr);

  va_end(argptr);

  return ier;
}

void PMMG2D_Init_parameters(PMMG2D_pParMesh parmesh, MPI_Comm comm) {
  int        flag;

  memset(&parmesh->info    ,0, sizeof(PMMG2D_Info));

  parmesh->info.mem                = PMMG2D_UNSET; /* [n/-1]   ,Set memory size to n Mbytes/keep the default value */
  parmesh->info.root               = PMMG2D_NUL;

  parmesh->ddebug                  = PMMG2D_NUL;
  parmesh->iter                    = PMMG2D_UNSET;
  parmesh->niter                   = PMMG2D_NITER;
  parmesh->info.fem                = MMG5_FEM;
  parmesh->info.ifc_layers         = PMMG2D_MVIFCS_NLAYERS;
  parmesh->info.nobalancing        = MMG5_OFF;
  parmesh->info.loadbalancing_mode = PMMG2D_LOADBALANCING_metis;
  parmesh->info.contiguous_mode    = PMMG2D_CONTIG_DEF;
  parmesh->info.sethmin            = PMMG2D_NUL;
  parmesh->info.sethmax            = PMMG2D_NUL;
  parmesh->info.fmtout             = PMMG2D_FMT_Unknown;
  parmesh->info.optim_interp       = PMMG2D_OPTIM_INTERP;
  parmesh->info.ratio_load_balance = PMMG2D_RATIO_LOAD_BALANCE;
  parmesh->info.interp_layers      = PMMG2D_INTERP_NLAYERS;

  /* Init MPI data */
  parmesh->comm   = comm;

  MPI_Initialized(&flag);
  parmesh->size_shm = 1;
  if ( flag ) {
    MPI_Comm_size( parmesh->comm, &parmesh->nprocs );
    MPI_Comm_rank( parmesh->comm, &parmesh->myrank );
  }
  else {
    parmesh->nprocs = 1;
    parmesh->myrank = PMMG2D_NUL;
  }

  /* ParMmg verbosity */
  if ( parmesh->myrank==parmesh->info.root ) {
    parmesh->info.imprim = PMMG2D_IMPRIM;
  }
  else {
    parmesh->info.imprim = PMMG2D_UNSET;
  }
  parmesh->info.imprim0  = PMMG2D_IMPRIM;

  /* Set the mmg verbosity to -1 (no output) */
  parmesh->info.mmg_imprim = PMMG2D_MMG_IMPRIM;
  parmesh->mesh->info.imprim = MG_MIN ( parmesh->info.imprim, PMMG2D_MMG_IMPRIM );

  // Set the nosizereq option for mmg. If nosizreq = 0, MMG runs into severe issues
  // due to the parallel interface boundaries in PMMG2D
  parmesh->mesh->info.nosizreq = PMMG2D_MMG_NOSIZREQ;

  // Set the hgradreq option for mmg. If hgradreq > 0, MMG runs into severe issues
  // due to the parallel interface boundaries in PMMG2D
  parmesh->mesh->info.hgradreq = PMMG2D_MMG_HGRADREQ;

  /* Default memory */
  PMMG2D_parmesh_SetMemGloMax( parmesh );
  PMMG2D_parmesh_SetMemMax( parmesh );

}

/**
 * \param parmesh pointer toward a parmesh structure for memory count.
 * \param buffer pointer toward the buffer to store file name.
 * \param name string to store into \a buffer.
 * \param defname default string for \a buffer if \a name is empty.
 * \return 1 if success, 0 if fail
 *
 * Set the file name \a name into the \a buffer string. If \a name is empty, use
 * \a defname as default file name. If \a defname is empty too, the buffer
 * remains non allocated.
 *
 * \warning for internal use only.
 */
int PMMG2D_Set_name(PMMG2D_pParMesh parmesh,char **buffer,
                    const char* name, const char* defname) {

  if ( *buffer ) {
    PMMG2D_DEL_MEM(parmesh,*buffer,char,"buffer unalloc");
  }

  if ( name && strlen(name) ) {
    PMMG2D_MALLOC(parmesh,*buffer,strlen(name)+1,char,"name",return 0);
    strcpy(*buffer,name);
  }
  else if ( defname ) {
    PMMG2D_MALLOC(parmesh,*buffer,strlen(defname)+1,char,"defname",return 0);
    strcpy(*buffer,defname);
  }
  /* Remark: for solution fields: a non allocated buffer allows to detect that
   * user don't provide input fields */

  return 1;
}

int PMMG2D_Set_iparameter(PMMG2D_pParMesh parmesh, int iparam,int val) {
  size_t      mem;

  switch ( iparam ) {
  case PMMG2D_IPARAM_verbose :
    if ( parmesh->myrank==parmesh->info.root ) {
      parmesh->info.imprim = val;
    }
    parmesh->info.imprim0 = val;

    break;
  case PMMG2D_IPARAM_mmgVerbose :
    parmesh->info.mmg_imprim = val;

    /* Set the mmg verbosity */
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_verbose,val) ) return 0;
    break;

  case PMMG2D_IPARAM_mem :
    if ( val <= 0 ) {
      fprintf( stdout,
        "  ## Warning: maximal memory per process must be strictly positive.\n");
      fprintf(stdout,"  Reset to default value.\n");
    } else {
      parmesh->info.mem = val;
    }
    PMMG2D_parmesh_SetMemGloMax(parmesh);
    parmesh->memMax = parmesh->memGloMax;
    mem = parmesh->memGloMax;

    /* Mesh reallocation if needed */
    if ( (parmesh->mesh->memCur >> MMG5_BITWIZE_MB_TO_B) > mem ) {
      fprintf(stderr,"\n  ## Error: %s: Maximal memory must be setted "
              "before reading the mesh.\n",__func__);
      return 0;
    }

    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_mem,(int)mem) ) return 0;
    
    break;
  case PMMG2D_IPARAM_nobalancing :
    parmesh->info.nobalancing = val;
    break;
  case PMMG2D_IPARAM_ifcLayers :
    parmesh->info.ifc_layers = val;
    break;
  case PMMG2D_IPARAM_interp_layers :
    parmesh->info.interp_layers = val;
    break;
  case PMMG2D_IPARAM_niter :
    parmesh->niter = val;
    break;
  case PMMG2D_IPARAM_debug :
    parmesh->ddebug = val;
    break;

  case PMMG2D_IPARAM_distributedOutput :

    if ( val == 2 ) {
      parmesh->info.fmtout = PMMG2D_FMT_Centralized_Distributed;
    }
    else if ( val == 1 ) {
      parmesh->info.fmtout = PMMG2D_FMT_Distributed;
    }
    else if ( val == 0 ) {
      parmesh->info.fmtout = PMMG2D_FMT_Centralized;
    }
    break;

  case PMMG2D_IPARAM_mmgDebug :

    /* Set the mmg debug mode */
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_debug,val) ) return 0;
    break;

  case PMMG2D_IPARAM_angle :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_angle,val) ) return 0;
    break;
  case PMMG2D_IPARAM_iso :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_iso,val) ) return 0;
    break;
  case PMMG2D_IPARAM_lag :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_lag,val) ) return 0;
    break;

  case PMMG2D_IPARAM_nofem :
    parmesh->info.fem    = (val==1)? 0 : 1;
    if ( !MMG2D_Set_iparameter(parmesh->mesh,parmesh->met,MMG2D_IPARAM_nofem,val) ) return 0;
    break;
  case PMMG2D_IPARAM_opnbdy :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_opnbdy,val) ) return 0;
    if( val ) {
      fprintf(stderr," ## Warning: Surface adaptation not supported with opnbdy."
          "\nSetting nosurf on.\n");
      if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_nosurf,val) ) return 0;
    }
    break;
  case PMMG2D_IPARAM_optim :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_optim,val) ) return 0;
    break;
  case PMMG2D_IPARAM_noinsert :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_noinsert,val) ) return 0;
    break;
  case PMMG2D_IPARAM_noswap :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_noswap,val) ) return 0;
    break;
  case PMMG2D_IPARAM_nomove :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_nomove,val) ) return 0;
    break;
  case PMMG2D_IPARAM_nosurf :
    if( !val && parmesh->mesh->info.opnbdy )
      fprintf(stderr," ## Warning: Surface adaptation not supported with opnbdy."
        "\nCannot set nosurf off.\n");
    else if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_nosurf,val) ) return 0;
    break;
  case PMMG2D_IPARAM_numberOfLocalParam :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_numberOfLocalParam,val) ) return 0;
    break;
  case PMMG2D_IPARAM_anisosize :
    if ( !MMG2D_Set_iparameter(parmesh->mesh,parmesh->met,MMG2D_IPARAM_anisosize,val) ) return 0;
    break;
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    return 0;
  }

  return 1;
}

int PMMG2D_Set_dparameter(PMMG2D_pParMesh parmesh, int dparam,double val){

  switch ( dparam ) {
  case PMMG2D_DPARAM_angleDetection :
    if ( !MMG2D_Set_dparameter(parmesh->mesh,NULL,MMG2D_DPARAM_angleDetection,val) ) return 0;
    break;
  case PMMG2D_DPARAM_hmin :
    parmesh->info.sethmin  = 1;
    if ( !MMG2D_Set_dparameter(parmesh->mesh,NULL,MMG2D_DPARAM_hmin,val) ) return 0;
    break;
  case PMMG2D_DPARAM_hmax :
    parmesh->info.sethmax  = 1;
    if ( !MMG2D_Set_dparameter(parmesh->mesh,NULL,MMG2D_DPARAM_hmax,val) ) return 0;
    break;
  case PMMG2D_DPARAM_hsiz :
    if ( !MMG2D_Set_dparameter(parmesh->mesh,NULL,MMG2D_DPARAM_hsiz,val) ) return 0;
    break;
  case PMMG2D_DPARAM_hgrad :
    if ( !MMG2D_Set_dparameter(parmesh->mesh,NULL,MMG2D_DPARAM_hgrad,val) ) return 0;
    break;
  case PMMG2D_DPARAM_hgradreq :
    if ( !MMG2D_Set_dparameter(parmesh->mesh,NULL,MMG2D_DPARAM_hgradreq,val) ) return 0;
    break;
  case PMMG2D_DPARAM_hausd :
    if ( val <=0 ) {
      fprintf(stderr,"\n  ## Error: %s: hausdorff number must be strictly"
              " positive.\n",__func__);
      return 0;
    }
    else {
        if ( !MMG2D_Set_dparameter(parmesh->mesh,NULL,MMG2D_DPARAM_hausd,val) ) return 0;
    }
    break;
  case PMMG2D_DPARAM_ls :
    if ( !MMG2D_Set_dparameter(parmesh->mesh,NULL,MMG2D_DPARAM_ls,val) ) return 0;
    break;
  case PMMG2D_DPARAM_load_balance :
    parmesh->info.ratio_load_balance = val;
    break;
  default :
    fprintf(stderr,"  ## Error: unknown type of parameter\n");
    return 0;
  }

  return 1;
}

/**
 * \param parmesh pointer toward a parmesh structure for memory count.
 * \param ghosts pointer toward the neighbouring partitions of each triangle.
 * \return 1 if success, 0 if fail
 *
 * Store the neighbour partitions of each triangle to identify ghost elements
 *
 */
int PMMG2D_Get_ghosts(PMMG2D_pParMesh parmesh, int** ghosts, int* size)
{
  MMG5_pMesh mesh;
  MMG5_pTria pt, pt_neigh;
  int ier = 1, ieresult;

  mesh = parmesh->mesh;

  // Store the interface nodes
  // Greatest reference number
  int max_ref = 0;
  for (int i = 1; i <= mesh->np; i++) {
    if (mesh->point[i].ref > max_ref) max_ref = mesh->point[i].ref;
  }

  int* interface_nodes = (int*) malloc((max_ref+1) * sizeof(int));
  for (int k = 0; k <= max_ref; k++) interface_nodes[k] = -1;
  for (int k = 0; k < parmesh->nip; k++)
    interface_nodes[mesh->point[(&parmesh->list_interface_nodes[k])->index[0]].ref] = k;

  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) )  continue;

    size[k-1] = 0;

    for (int i = 0; i < 3; i++) {

      int k1 = interface_nodes[mesh->point[pt->v[i]].ref];

      if (k1 == -1) continue;

      for (int j = 1; j < (&parmesh->list_interface_nodes[k1])->nb; j++) {

        int should_add = 1;

        for (int l = 0; l < size[k-1]; l++) {
          if ((&parmesh->list_interface_nodes[k1])->proc[j] == ghosts[k-1][l]) {
            should_add = 0;
            break;
          }
        }

        if (should_add) {
          size[k-1]++;
          ghosts[k-1] = realloc(ghosts[k-1], size[k-1] * sizeof(int));
          ghosts[k-1][size[k-1]-1] = (&parmesh->list_interface_nodes[k1])->proc[j];
        }

      }

    }

  }

  free(interface_nodes);

  return 1;
}
