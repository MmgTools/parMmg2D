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
 * \file zaldy_pmmg2d.c
 * \brief Memory management
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Memory management for PARMMG2D library.
 *
 */

#include "parmmg2d.h"

/**
 * \param parmesh pointer to pmmg structure
 * \param memReq  size of memory in Mb. If memReq is zero then it is
 *                automatically set to half of the machine's available memory.
 *                On machines with multicore processors
 *                the total available memory is shared equally to pmmg processes
 *                running on the same machine.
 *                If memReq is negative or more than the detected available
 *                memory, then the requested value is discarded and the maximum
 *                allowed memory is set to half the detected available memory.
 *
 *  Sets the maximum amount of memory that a parmmg process is allowed
 *  to use depending on the user specifications and the available
 *  memory. On machines with multicore processors the memory is shared
 *  equally to pmmg processes running on the same machine. 
 */
void PMMG2D_parmesh_SetMemGloMax( PMMG2D_pParMesh parmesh )
{
  size_t   maxAvail = 0;
  MPI_Comm comm_shm = 0;
  int      flag;

  assert ( (parmesh != NULL) && "trying to set glo max mem in empty parmesh" );

  /** Step 1: Get the numper of processes per node */
  MPI_Initialized( &flag );

  if ( flag ) {
    MPI_Comm_split_type( parmesh->comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                         &comm_shm );
    MPI_Comm_size( comm_shm, &parmesh->size_shm );
  }
  else {
    parmesh->size_shm = 1;
  }

  /** Step 2: Set maximal memory per process depending on the -m option setting */
  maxAvail = MMG5_memSize();

  if ( parmesh->info.mem <= 0 ) {
    /* Nos users specifications */
    if ( !maxAvail ) {
      /* default value when not able to compute the available memory = 800 MB */
      printf("  Maximum memory per process set to default value: %d MB.\n",MMG5_MEMMAX);
      parmesh->memGloMax = MMG5_MEMMAX << 20;
    }
    else {
      /* maximal memory = total physical memory */
      parmesh->memGloMax = maxAvail;
    }
  }
  else {
    /* memory asked by user if possible, otherwise total physical memory */
    if ( maxAvail && (size_t)parmesh->info.mem*MMG5_MILLION > maxAvail ) {
      fprintf(stderr,"\n  ## Warning: %s: asking for %d MB of memory per process ",
              __func__,parmesh->info.mem);
      fprintf(stderr,"when only %zu available.\n",maxAvail/MMG5_MILLION);
    }
    else {
      parmesh->memGloMax= (size_t)parmesh->info.mem*MMG5_MILLION;
    }
  }

  if ( abs(parmesh->info.imprim) > 4 || parmesh->ddebug ) {
    fprintf(stdout,"  MAXIMUM MEMORY AUTHORIZED PER PROCESS (MB)    %zu\n",
            parmesh->memGloMax/MMG5_MILLION);
  }
}

/**
 * \param parmesh parmesh structure
 *
 * \return 1 if success, 0 if fail
 *
 * Set the maximum memory that parmesh can use.
 */
int PMMG2D_parmesh_SetMemMax( PMMG2D_pParMesh parmesh ) {

  parmesh->memMax = parmesh->memGloMax;
  parmesh->mesh->memMax = parmesh->memGloMax;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param met pointer toward the metric structure
 *
 * \return 0 if fail, 1 otherwise
 *
 * memory repartition between the point, xpoint, tetra, xtetra arrays for the
 * memMax amout of memory available.
 *
 */
static inline
int PMMG2D_memOption_memRepartition(MMG5_pMesh mesh,MMG5_pSol met) {
  size_t     usedMem,avMem;
  size_t     npadd;
  int        ctri,bytes;

  // Compute the needed initial memory
  usedMem = (mesh->np+1)*sizeof(MMG5_Point)
    + (mesh->xp+1)*sizeof(MMG5_xPoint) + (mesh->nt+1)*sizeof(MMG5_Tria) + (mesh->na+1)*sizeof(MMG5_pEdge);

  if ( mesh->adja )
    usedMem += (3*mesh->nt+1)*sizeof(int);

  if ( met->m )
    usedMem += met->size*(mesh->np+1)*sizeof(double);

  if ( usedMem > mesh->memMax  ) {
    fprintf(stderr,"\n  ## Error: %s: %zu Mo of memory ",__func__,mesh->memMax/MMG5_MILLION);
    fprintf(stderr,"is not enough to load mesh. You need to ask %zu Mo per process minimum.\n",
            usedMem/MMG5_MILLION+1);
    fprintf(stderr,"\nTry to use the -m option to impose the maximal memory per process.\n");
    return 0;
  }


  /** Try to estimate the memody usage of adding a point (with the related
   *  tria, solution...) in order to reduce the maximum size when
   *  possible. */

  ctri = 2;
  /* Euler-poincare: ne = 6*np; nt = 2*np; na = np/5 *
   * point+tria+tets+adja+adjt+sol+item */
  bytes = sizeof(MMG5_Point) + sizeof(MMG5_xPoint) +
          ctri*sizeof(MMG5_Tria) ;

  if ( mesh->adja )
    bytes += 4*6*sizeof(int);

  if ( met->m )
    bytes += met->size*sizeof(double);

  avMem = mesh->memMax-usedMem;

  /* The number of points that can be added is approximately given by the
   * ratio between the available memory and the memory usage of a
   * point+related structures */
  npadd = (size_t) ( (double)avMem/bytes );

  /* Shrink size if too big */
  mesh->npmax = MG_MIN(mesh->npmax,mesh->np+npadd);
  mesh->xpmax = MG_MIN(mesh->xpmax,mesh->xp+npadd);
  mesh->ntmax = MG_MIN(mesh->ntmax,ctri*npadd+mesh->nt);

  met->npmax  = mesh->npmax;

  if ( abs(mesh->info.imprim) > 5 || mesh->info.ddebug ) {
    fprintf(stdout,"  MAXIMUM MEMORY PER PROCESS AUTHORIZED (Mo)    %zu\n",
            mesh->memMax/MMG5_MILLION);

    fprintf(stdout,"  MMG2D_NPMAX    %d\n",mesh->npmax);
    fprintf(stdout,"  MMG2D_XPMAX    %d\n",mesh->xpmax);
    fprintf(stdout,"  MMG2D_NTMAX    %d\n",mesh->ntmax);
  }

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Keep track of empty links for the triangles and points array.
 *
 */
int PMMG2D_link_mesh( MMG5_pMesh mesh ) {
  int k;

  /* keep track of empty links */
  if ( mesh->npmax > mesh->np ) {
    mesh->npnil = mesh->np + 1;
    for (k=mesh->npnil; k<mesh->npmax-1; k++)
      mesh->point[k].tmp  = k+1;
  }
  else {
    assert ( mesh->np == mesh->npmax );
    mesh->npnil = 0;
  }

  if ( mesh->ntmax > mesh->nt ) {
    mesh->nenil = mesh->nt + 1;
    for (k=mesh->nenil; k<mesh->ntmax-1; k++)
      mesh->tria[k].v[2] = k+1;
  }
  else {
    assert ( mesh->nt == mesh->ntmax );
    mesh->nenil = 0;
  }
  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Allocation of the array fields of the mesh for the given npmax, xpmax, nemax,
 * xtmax.
 *
 */
int PMMG2D_setMeshSize_alloc( MMG5_pMesh mesh ) {

  PMMG2D_CALLOC(mesh,mesh->point,mesh->npmax+1,MMG5_Point,
                "vertices array", return 0);

  PMMG2D_CALLOC(mesh,mesh->xpoint,mesh->xpmax+1,MMG5_xPoint,
                "boundary vertices array", return 0);

  if ( mesh->nt ) {
    PMMG2D_CALLOC(mesh,mesh->tria,mesh->nt+1,MMG5_Tria,
                  "triangles array", return 0);
  }

  return ( PMMG2D_link_mesh( mesh ) );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param npmax_old old maximum number of points.
 * \param xpmax_old old maximum number of boundary points.
 * \param ntmax_old old maximum number of triangles.
 * \param namax_old old maximum number of edges.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Reallocation of the array fields of the mesh for the given npmax,
 * xpmax, ntmax, xtmax.
 *
 */
int PMMG2D_setMeshSize_realloc( MMG5_pMesh mesh,int npmax_old,int xpmax_old, int ntmax_old, int namax_old ) {

  PMMG2D_RECALLOC(mesh,mesh->point,mesh->npmax+1,npmax_old+1,MMG5_Point,
                  "vertices array", return 0);

  PMMG2D_RECALLOC(mesh,mesh->xpoint,mesh->xpmax+1,xpmax_old+1,MMG5_xPoint,
                  "boundary vertices array", return 0);

  PMMG2D_RECALLOC(mesh,mesh->tria,mesh->ntmax+1,ntmax_old+1,MMG5_Tria,
                  "triangle array", return 0);

  PMMG2D_RECALLOC(mesh,mesh->edge,mesh->namax+1,namax_old+1,MMG5_Edge,
                  "edge array", return 0);

  if ( mesh->adja ) {
    PMMG2D_RECALLOC(mesh,mesh->adja,3*mesh->ntmax+5,3*ntmax_old+5,int,
                    "adja array", return 0);
  }

  return ( PMMG2D_link_mesh( mesh ) );
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nt number of triangles.
 * \param xp number of boundary point
 * \param xt number of boundary tetra
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Check the input mesh size and assign their values to the mesh.
 *
 */
int PMMG2D_setMeshSize_initData(MMG5_pMesh mesh, int np, int nt, int xp, int xt ) 
{
  if ( ( (mesh->info.imprim > PMMG2D_VERB_DETQUAL) || mesh->info.ddebug ) &&
       ( mesh->point || mesh->xpoint ) )
    fprintf(stderr,"\n  ## Warning: %s: old mesh deletion.\n",__func__);

  if ( !np ) {
    fprintf(stderr,"  ** MISSING DATA:\n");
    fprintf(stderr,"     Your mesh must contains at least points.\n");
    return 0;
  }
  if ( !nt && (mesh->info.imprim > PMMG2D_VERB_DETQUAL || mesh->info.ddebug) ) {
    fprintf(stderr,"  ** WARNING:\n");
    fprintf(stderr,"     Your mesh don't contains tetrahedra.\n");
  }

  if ( mesh->point )
    MMG5_DEL_MEM(mesh,mesh->point);
  if ( mesh->tria )
    MMG5_DEL_MEM(mesh,mesh->tria);
  if ( mesh->quadra )
    MMG5_DEL_MEM(mesh,mesh->quadra);
  if ( mesh->edge )
    MMG5_DEL_MEM(mesh,mesh->edge);

  mesh->np  = np;
  mesh->nt  = nt;
  mesh->xp  = xp;
  mesh->xt  = xt;

  mesh->npi = mesh->np;
  mesh->nti = mesh->nt;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param ne number of tetrahedra.
 * \param nt number of triangles.
 * \param xp number of boundary points.
 * \param xt number of boundary tetra.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Check the input mesh size and assign their values to the mesh.
 *
 */
int PMMG2D_setMeshSize(MMG5_pMesh mesh, int np, int nt, int xp, int xt ) {

  /* Check input data and set mesh->np/nt to the suitable values */
  if ( !PMMG2D_setMeshSize_initData(mesh,np,nt,xp,xt) )
    return 0;

  mesh->npmax  = mesh->np;
  mesh->ntmax  = mesh->nt;
  mesh->xpmax  = mesh->xp;
  mesh->xtmax  = mesh->xt;

  /* Mesh allocation and linkage */
  if ( !PMMG2D_setMeshSize_alloc( mesh ) ) return 0;

  return 1;
}


/**
 * \param parmesh parmesh structure to adjust
 * \param fitMesh if 1, set maximum mesh size at its exact size.
 *
 * \return 1 if success, 0 if fail
 *
 * Update the size of the mesh.
 *
 */
int PMMG2D_updateMeshSize( PMMG2D_pParMesh parmesh, int fitMesh )
{
  MMG5_pMesh mesh;
  MMG5_pSol  met,ls,disp;
  int        npmax_old,xpmax_old,ntmax_old,namax_old;

  mesh = parmesh->mesh;

  /* Force the MMG5_memOption_memSet function to find the wanted memMax value
   * in MMG2D_Set_meshSize, MMG2D_Set_iparameter, MMG2D_zaldy (for i/o) */
  mesh->info.mem = parmesh->memGloMax/MMG5_MILLION;

  /* Memory repartition for the MMG meshes arrays */
  npmax_old = mesh->npmax;
  xpmax_old = mesh->xpmax;
  ntmax_old = mesh->ntmax;
  namax_old = mesh->namax;
  if ( fitMesh ) {
    mesh->npmax = mesh->np;
    mesh->xpmax = mesh->xp;
    mesh->ntmax = mesh->nt;
    mesh->namax = mesh->na;
  }
  else {
    mesh->npmax = 1.5*mesh->np;
    mesh->xpmax = 1.5*mesh->xp;
    mesh->ntmax = 1.5*mesh->nt;
    mesh->namax = 1.5*mesh->na;
  }

  met = parmesh->met;
  if ( !PMMG2D_memOption_memRepartition(mesh,met) ) return 0;

  if ( !PMMG2D_setMeshSize_realloc(mesh,npmax_old,xpmax_old,ntmax_old,namax_old) )
    return 0;

  if ( met ) {
    met->np    = mesh->np;
    met->npmax = mesh->npmax;
  }
  ls = parmesh->ls;
  if ( ls ) {
    ls->np      = mesh->np;
    ls->npmax   = mesh->npmax;
  }

  disp = parmesh->disp;
  if ( disp ) {
    disp->np    = mesh->np;
    disp->npmax = mesh->npmax;
  }

  if ( met && met->m )
    PMMG2D_RECALLOC(mesh,met->m,met->size*(met->npmax+1),met->size*(npmax_old+1),
                   double,"metric array",return 0);

  if ( ls && ls->m ) {
    PMMG2D_RECALLOC(mesh,ls->m,ls->size*(ls->npmax+1),
                   ls->size*(npmax_old+1),double,"ls_array",
                   return 0);
  }

  if ( disp && disp->m ) {
    PMMG2D_RECALLOC(mesh,disp->m,disp->size*(disp->npmax+1),
                   disp->size*(npmax_old+1),double,"disp_array",
                   return 0);
  }

  return 1;
}

