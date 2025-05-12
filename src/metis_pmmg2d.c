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
 * \file metis_pmmg2d.c
 * \brief Partition mesh using metis
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */
#include "metis_pmmg2d.h"
#include "linkedlist_pmmg2d.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param xadj array of shifts (CSR).
 * \param adjncy array of adjacents (CSR).
 * \param adjwgt array of weights (CSR).
 * \param filename filename prefix.
 *
 * \return 1 if success, 0 if fail.
 *
 * Save mesh graph and weight on file in Medit format.
 *
 */
int PMMG2D_saveGraph( PMMG2D_pParMesh parmesh, const char *filename, idx_t *part)
{
  MMG5_pMesh mesh = parmesh->mesh;
  MMG5_pTria pt;
  MMG5_pPoint ppt;
  char *sname,*smesh;
  FILE *fmesh;
  int k;

  PMMG2D_CALLOC(parmesh,sname,strlen(filename)+9,char,"file name prefix",return 0);
  PMMG2D_CALLOC(parmesh,smesh,strlen(filename)+15,char,"mesh file name",return 0);
  sprintf(sname,"%s-P%02d-I%02d",filename,parmesh->myrank,parmesh->iter);
  strcpy(smesh,sname);
  strcat(smesh,".mesh");

  /* Open files and write headers */
  fmesh = fopen(smesh,"w");
  fprintf(fmesh,"MeshVersionFormatted 2\n");
  fprintf(fmesh,"\nDimension 2\n");

  /* Write vertices */
  fprintf(fmesh,"\nVertices\n%d\n",mesh->np);
  for( k = 1; k <= mesh->np; k++ ) {
    ppt = &mesh->point[k];
    fprintf(fmesh,"%f %f %d\n",
            ppt->c[0],
            ppt->c[1],
            ppt->ref);
  }

  /* Write triangles and solution on triangles */
  fprintf(fmesh,"\nTriangles\n%d\n",mesh->nt);
  for( k = 1; k <= mesh->nt; k++ ) {
    pt   = &mesh->tria[k];
    fprintf(fmesh,"%d %d %d %d\n",
            pt->v[0],
            pt->v[1],
            pt->v[2],
            part[k-1]);
  }

  /* Close files */
  fprintf(fmesh,"\n\nEnd");
  fclose(fmesh);

  /* Free memory */
  PMMG2D_DEL_MEM(parmesh,sname,char,"file name prefix");
  PMMG2D_DEL_MEM(parmesh,smesh,char,"mesh file name");
  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param part    elements partition array
 * \param ne      nb of elements
 * \param nproc   nb of groups for partitioning
 *
 * \return 1 if no empty partitions or successfully corrected, 0 if fail
 *
 * Check if metis has returned empty partitions, correct partitioning if so.
 *
 */
int PMMG2D_correct_meshElts2metis( PMMG2D_pParMesh parmesh,idx_t* part,idx_t ne,idx_t nproc ) {
  PMMG2D_valLnkdList **partlist;
  idx_t            iproc,ie;
  int              nempt,iempt;


  /* Initialize lists */
  PMMG2D_CALLOC(parmesh, partlist, nproc, PMMG2D_valLnkdList*,"array of list pointers",return 0);
  for( iproc=0; iproc<nproc; iproc++ ) {
    PMMG2D_CALLOC(parmesh,partlist[iproc],1,PMMG2D_valLnkdList,"linked list pointer",return 0);
    if( !PMMG2D_valLnkdListNew(parmesh,partlist[iproc],iproc,PMMG2D_LISTSIZE) ) return 0;
  }

  /* Fill the lists */
  assert ( ne >= nproc &&  "not enough elements for the number of partitions" );

  for( ie = 0; ie < ne; ie++ ) {
    iproc = part[ie];
    if( !PMMG2D_add_val2lnkdList(parmesh,partlist[iproc],ie) ) return 0;
  }

  /* Sort lists based on nb. of entities, in ascending order */
  qsort(partlist,nproc,sizeof(PMMG2D_valLnkdList*),PMMG2D_compare_valLnkdListLen);

  /* Count empty partitions */
  nempt = 0;
  for( iproc = 0; iproc < nproc; iproc++ ) {
    if( partlist[iproc]->nitem ) break;
    nempt++;
  }

  assert( nempt < nproc );
  if( !nempt ) {
    /* Deallocate lists and return */
    for( iproc = 0; iproc < nproc; iproc++ ) {
      PMMG2D_DEL_MEM(parmesh,partlist[iproc]->item,PMMG2D_lnkdVal,"linked list array");
      PMMG2D_DEL_MEM(parmesh,partlist[iproc],PMMG2D_valLnkdList,"linked list pointer");
    }
    PMMG2D_DEL_MEM(parmesh,partlist,PMMG2D_valLnkdList*,"array of linked lists");

    return 1;
  }

  fprintf(stdout,"   ### Warning: Empty partitions on proc %d, nelts %d\n",parmesh->myrank,ne);
  /** Correct partitioning */
  assert ( nproc > 1 );
  iproc = nproc-1;
  iempt = 0;
  while( nempt ) {
    /* Get next "reservoir" proc */
    if ( iproc == nproc-1 ) {
      while( partlist[iproc]->nitem <= partlist[iproc-1]->nitem ) {

        /* list are sorted depending to their number of items so iproc has more
         * items than iproc-1 */
        assert ( partlist[iproc]->nitem == partlist[iproc-1]->nitem );
      iproc--;
      }
      iproc--;
    }
    ++iproc;

    /* if ne > nproc, normally, we can fill the empty procs without emptying a
     * proc with only 1 item */
    assert ( partlist[iproc]->nitem > 1 && "not enough elements for the"
             " number of partitions");

    /* Pop entity ie from iproc, add to iempt */
    if( !PMMG2D_pop_val_lnkdList(parmesh,partlist[iproc],&ie) ) return 0;
    if( !PMMG2D_add_val2lnkdList(parmesh,partlist[iempt],ie) ) return 0;
    /* Update partition table and go on to next empty proc */
    part[ie] = partlist[iempt]->id;
    iempt++;
    nempt--;
  }

  /* Deallocate lists */
  for( iproc=0; iproc<nproc; iproc++ ) {
    PMMG2D_DEL_MEM(parmesh,partlist[iproc]->item,PMMG2D_lnkdVal,"linked list array");
    PMMG2D_DEL_MEM(parmesh,partlist[iproc],PMMG2D_valLnkdList,"linked list pointer");
  }
  PMMG2D_DEL_MEM(parmesh,partlist,PMMG2D_valLnkdList*,"array of linked lists");

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param part pointer of an array containing the partitions (at the end)
 * \param nproc number of partitions asked
 *
 * \return  1 if success, 0 if fail
 *
 * Use metis to partition the mesh into nproc submeshes
 *
 */
int PMMG2D_part_meshElts2metis( PMMG2D_pParMesh parmesh, idx_t* part, idx_t nproc )
{
  idx_t      *xadj,*adjncy,*vwgt,*adjwgt;
  idx_t      adjsize;
  idx_t      nelt = parmesh->mesh->nt;
  idx_t      ncon = 1; // number of balancing constraint
  idx_t      options[METIS_NOPTIONS];
  idx_t      objval = 0;
  int        ier = 0;
  int        status = 1;

  xadj = adjncy = vwgt = adjwgt = NULL;

  /* Set contiguity of partitions if using Metis also for graph partitioning */
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_CONTIG] = ( parmesh->info.contiguous_mode &&
    (parmesh->info.loadbalancing_mode & PMMG2D_LOADBALANCING_metis) );

  /** Build the graph */
  if ( !PMMG2D_graph_meshElts2metis(parmesh,parmesh->mesh,parmesh->met,&xadj,&adjncy,&adjwgt,&adjsize) )
    return 0;

  /** Call metis and get the partition array */
  if( nproc >= 8 ) {
    ier = METIS_PartGraphKway( &nelt,&ncon,xadj,adjncy,vwgt,NULL,adjwgt,&nproc,
                               NULL,NULL,options,&objval, part );
  }
  else {
    ier = METIS_PartGraphRecursive( &nelt,&ncon,xadj,adjncy,vwgt,NULL,adjwgt,&nproc,
                               NULL,NULL,options,&objval, part );
  }

  if ( ier != METIS_OK ) {
    switch ( ier ) {
      case METIS_ERROR_INPUT:
        fprintf(stderr, "METIS_ERROR_INPUT: input data error\n" );
        break;
      case METIS_ERROR_MEMORY:
        fprintf(stderr, "METIS_ERROR_MEMORY: could not allocate memory error\n" );
        break;
      case METIS_ERROR:
        fprintf(stderr, "METIS_ERROR: generic error\n" );
        break;
      default:
        fprintf(stderr, "METIS_ERROR: update your METIS error handling\n" );
        break;
    }
    status = 0;
  }

  /** Correct partitioning to avoid empty partitions */
  if( !PMMG2D_correct_meshElts2metis( parmesh,part,nelt,nproc ) ) return 0;

  PMMG2D_DEL_MEM(parmesh, adjwgt, idx_t, "deallocate adjwgt" );
  PMMG2D_DEL_MEM(parmesh, adjncy, idx_t, "deallocate adjncy" );
  PMMG2D_DEL_MEM(parmesh, xadj, idx_t, "deallocate xadj" );
  PMMG2D_DEL_MEM(parmesh, vwgt, idx_t, "deallocate vwgt" );

  return status;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param met  pointer toward the met structure.
 * \param pt1  pointer toward the first triangle.
 * \param pt2  pointer toward the second triangle.
 *
 * \return The weight value
 *
 * Compute an element weight to be used for the metis weight on parallel
 * interfaces.
 *
 */
double PMMG2D_computeWgt( MMG5_pMesh mesh, MMG5_pSol met, MMG5_pTria pt1, MMG5_pTria pt2 ) {
  double       len, res, alpha;
  double       m1[3],m2[3];
  int          i,j,ip1,ip2;
  MMG5_pPoint  pp1,pp2;

  ip1 = -1;
  ip2 = -1;
  alpha = 28.;

  if ( met && met->m ) {
    // First find the two common vertices between the two triangles
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        if (pt1->v[i] == pt2->v[j]) {
          if (ip1 == -1) {
            ip1 = pt1->v[i];
            break;
          }
          else {
            ip2 = pt1->v[i];
            break;
          }
        }
      }
    }

    if ( ip2 == -1 || ip1 > mesh->np || ip2 > mesh->np ) {
      fprintf( stderr,"  ## Error in PMMG2D_computeWgt: likely due to bad triangle adjacency.\n" );
      return 1.;
    }


    pp1 = &mesh->point[ip1];
    pp2 = &mesh->point[ip2];

    // Give the metrics at both points and their coordinates to compute the edge length
    for ( i = 0; i < 3; ++i ) m1[i] = met->m[3*ip1+i];
    for ( i = 0; i < 3; ++i ) m2[i] = met->m[3*ip2+i];

    // Compute edge length
    len = long_ani(pp1->c, pp2->c, m1, m2);
    if( len <= 1.0 ) {
      res = len-1.0;
    }
    else {
      res = 1.0/len-1.0;
    }

    res = MG_MIN(1.0/exp(alpha*res),PMMG2D_WGTVAL_HUGEINT);
  }
  else {
    res = PMMG2D_WGTVAL_HUGEINT;
  }

  return res;
}

/**
 * \param parmesh pointer toward the PMMG parmesh structure
 * \param mesh pointer toward a MMG5 mesh structure
 * \param xadj pointer toward the position of the elt adjacents in adjncy
 * \param adjncy pointer toward the list of the adjacent of each elt
 * \param nadjncy number of data in adjncy array
 *
 * \return  1 if success, 0 if fail
 *
 * Build the metis graph with the mesh elements as metis nodes.
 *
 * \warning the mesh must be packed
 *
 */
int PMMG2D_graph_meshElts2metis( PMMG2D_pParMesh parmesh,MMG5_pMesh mesh,MMG5_pSol met,
                                 idx_t **xadj,idx_t **adjncy,idx_t **adjwgt,
                                 idx_t *nadjncy ) {
  MMG5_pTria   pt;
  int          *adjt;
  int          j,k,iadr,jel,count,nbAdj,wgt,ier;

  // Step 1: mesh adjacency creation
  if ( (!mesh->adja) && (1 != MMG2D_hashTria( mesh )) ) {
    fprintf( stderr,"  ## PMMG2D Hashing problem (1).\n" );
    return 0;
  }

  // Step 2: build the metis graph 
  PMMG2D_CALLOC(parmesh, (*xadj), mesh->nt+1, idx_t, "allocate xadj",
                return 0);

  // 1) Count the number of adjacent triangles of each elements and fill xadj
  (*xadj)[0] = 0;
  (*nadjncy) = 0;
  for ( k = 1; k <= mesh->nt; k++ ) {
    nbAdj = 0;
    iadr = 3*(k-1) + 1;
    adjt = &mesh->adja[iadr];
    for ( j = 0; j < 3; j++ )
      if ( adjt[j] ) nbAdj++;

    (*nadjncy) += nbAdj;
    (*xadj)[k] = (*nadjncy);
  }

  // 2) List the adjacent of each elts in adjncy
  ier = 1;
  ++(*nadjncy);
  PMMG2D_CALLOC(parmesh, (*adjncy), (*nadjncy), idx_t, "allocate adjncy", ier=0;);
  if( !ier ) {
    PMMG2D_DEL_MEM(parmesh, (*xadj), idx_t, "deallocate xadj" );
    return ier;
  }

  // Don't compute weights at mesh distribution, or if output load balancing is required at last iter
  if( (parmesh->iter != PMMG2D_UNSET) &&
      ((parmesh->iter < parmesh->niter-1) || parmesh->info.nobalancing) ) {
    PMMG2D_CALLOC(parmesh, (*adjwgt), (*nadjncy), idx_t, "allocate adjwgt", ier=0;);
    if( !ier ) {
      PMMG2D_DEL_MEM(parmesh, (*xadj), idx_t, "deallocate xadj" );
      PMMG2D_DEL_MEM(parmesh, (*adjncy), idx_t, "deallocate adjncy" );
      return ier;
    }
  }

  count = 0;
  for( k = 1; k <= mesh->nt; k++ ) {
    iadr = 3*(k-1) + 1;
    adjt = &mesh->adja[iadr];
    pt   = &mesh->tria[k];
    for ( j = 0; j < 3; j++ ) {
      jel = adjt[j] / 3;
      if ( !jel ) continue;

      // Assign graph edge weights
      if( *adjwgt ) {
        // Compute weight using edge size
        wgt = (int)PMMG2D_computeWgt(mesh,met,pt,&mesh->tria[jel]);

        (*adjwgt)[count] = MG_MAX(wgt,1);
      }

      (*adjncy)[count++]   = jel-1;
    }
    assert( count == ( (*xadj)[k] ) );
  }

  return ier;
}

/**
 * \param parmesh pointer toward a PMMG parmesh structure.
 *
 * \return 0 (on all procs) if fail, 1 otherwise
 *
 * Send the local meshes to the corresponding procs.
 */
int PMMG2D_Metis_exchange(PMMG2D_pParMesh parmesh)
{
  int ier, ieresult;
  int npt_global, npp_global;

  MMG5_pTria pt_global = NULL;
  MMG5_pPoint pp_global = NULL;
  double* mm_global = NULL;
  MMG5_pTria pt;
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;

  int* list_deleted_triangles = (int*) malloc((mesh->nt+1) * sizeof(int));

  for (int i = 0; i < mesh->nt; i++) list_deleted_triangles[i] = 0;

  // 1) Find the triangles to exchange
  int n = 0;
  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;
    if (pt->base == parmesh->myrank) continue;
    list_deleted_triangles[n++] = k;
  }

  // 2) Exchange the triangles, points and metrics
  ier = PMMG2D_exchange_from_root( parmesh, mesh->nt, &pt_global, &pp_global, &mm_global, &npt_global, &npp_global );

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error in the mpi exchanges during the Metis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 3) Remove the sent triangles from the local mesh
  int j = 0;
  for (int k = 1; k <= mesh->nt; k++ ) {
    if (list_deleted_triangles[j] == k) {
      j += 1;
      continue;
    }
    mesh->tria[k-j] = mesh->tria[k];
  }

  free(list_deleted_triangles);

  mesh->nt -= j;
  int nt_tmp = mesh->nt - j; // temporary number of triangles, used to check the contiguity (removing twice the triangles)

  PMMG2D_REALLOC(mesh, mesh->tria, mesh->nt+1, mesh->ntmax+1,
                 MMG5_Tria,"Metis decomposition realloc triangles",return 0);
  PMMG2D_REALLOC(mesh, mesh->point, mesh->np+1, mesh->npmax+1,
                 MMG5_Point,"Metis decomposition realloc points", return 0);
  PMMG2D_REALLOC(mesh, parmesh->met->m, parmesh->met->size*(mesh->np+1), parmesh->met->size*(mesh->npmax+1),
                 double, "Metis decomposition realloc metrics", return 0);

  // 4) Remove the points that do not belong to the mesh anymore
  ier = PMMG2D_remove_points( parmesh );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when removing the points during the Metis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 5) Add the triangles to their corresponding new process
  ier = PMMG2D_add_triangles( parmesh, pt_global, pp_global, mm_global, npt_global, npp_global);
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when adding the new triangles during the Metis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  PMMG2D_update( parmesh );

  // 6) Untag the interface parallel points and edges
  PMMG2D_untag_parallel( parmesh );

  // 7) Mark the interface nodes
  ier = PMMG2D_mark_parallel_interface_nodes( parmesh );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when marking the interface nodes during the contiguity check in the Metis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Free the pointers
  PMMG2D_DEL_MEM(mesh, pt_global, MMG5_Tria, "Delete pt_global in the Metis decomposition");
  PMMG2D_DEL_MEM(mesh, pp_global, MMG5_Point, "Delete pp_global in the Metis decomposition");
  PMMG2D_DEL_MEM(mesh, mm_global, double, "Delete mm_global in the Metis decomposition");

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if it fails, 1 otherwise. 
 *
 * Decompose the mesh sequentially using Metis and spread it over processors.
 *
 */
int PMMG2D_Metis_decomposition(PMMG2D_pParMesh parmesh)
{
  int ier;
  idx_t      *part;
  MMG5_pMesh mesh;
  MMG5_pTria pt;
  mesh = parmesh->mesh;

  if (!parmesh->myrank && mesh->edge) {
    PMMG2D_DEL_MEM(mesh, mesh->edge, MMG5_Edge, "Delete edges when decomposing the mesh using Metis");
    PMMG2D_REALLOC(mesh, mesh->edge, 1, 0, MMG5_Edge, "Realloc edge when decomposing the mesh using Metis (1)", return 0);
    mesh->na = 0;
  }

  // There is nothing to distribute on just 1 proc
  if( parmesh->nprocs == 1 ) return 1;

  // Partitions the mesh on root processor
  if( parmesh->myrank == parmesh->info.root ) {

    // Attribute to ref the point index
    for ( int i = 1; i <= parmesh->mesh->np; i++) parmesh->mesh->point[i].ref = i;

#ifdef USE_METIS
    // Call metis for partionning
    PMMG2D_CALLOC ( parmesh,part,parmesh->mesh->nt,idx_t,"allocate metis buffer", ier=5 );

    // Call metis, or recover a custom partitioning if provided (only to debug
    // the interface displacement, adaptation will be blocked) 
    if( !PMMG2D_PREDEF_PART ) {
      if ( PMMG2D_part_meshElts2metis( parmesh, part, parmesh->nprocs ) ) {
        // Store the partition number of each triangle in base
        for (int k = 1; k <= parmesh->mesh->nt; k++)
          parmesh->mesh->tria[k].base = part[k-1];
      }
      else ier = 5;
    } else {
      for(int k = 1; k <= parmesh->mesh->nt; k++ ) {
        parmesh->mesh->tria[k].base = parmesh->mesh->tria[k].ref;
      }
    }

    PMMG2D_DEL_MEM(parmesh,part,idx_t,"deallocate metis buffer");
#endif
  }

  // Send the local meshes from root process to the corresponding procs
  ier = PMMG2D_Metis_exchange(parmesh);

  int ieresult;
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

  // Rebuild the boundary edges
  for (int k = 1; k <= mesh->nt; k++ ) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) continue;

    for (int i = 0; i < 3; i++ ) {

      if (pt->tag[i] & MG_BDY) {
        mesh->na++;
        PMMG2D_REALLOC(mesh, mesh->edge, mesh->na+1, mesh->na,
                       MMG5_Edge, "Realloc edge when decomposing the mesh using Metis (2)", return 0);
        int i1 = MMG5_inxt2[i];
        int i2 = MMG5_inxt2[i1]; 
        mesh->edge[mesh->na].a = pt->v[i1];
        mesh->edge[mesh->na].b = pt->v[i2];
        mesh->edge[mesh->na].tag = pt->tag[i];
        mesh->edge[mesh->na].ref = mesh->na;
        pt->edg[i] = mesh->na;
      }
    }
  }

  return ieresult;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param partTab list of partitions of each triangle
 * \return 0 if it fails, 1 otherwise. 
 *
 * Modify the local mesh to fit with the new domain decomposition
 *
 */
int PMMG2D_ParMetis_exchange( PMMG2D_pParMesh parmesh, int *partTab )
{
  int list_index_pt[parmesh->nprocs];
  int list_index_pp[parmesh->nprocs];

  MMG5_pTria pt_global = NULL;
  MMG5_pPoint pp_global = NULL;
  double* mm_global = NULL;
  MMG5_pTria pt;
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;

  MMG5_pTria* list_pt = (MMG5_pTria*) malloc(parmesh->nprocs * sizeof(MMG5_pTria));
  for (int k = 0; k < parmesh->nprocs; k++) list_pt[k] = (MMG5_pTria) malloc(mesh->nt * sizeof(MMG5_Tria));
  MMG5_pPoint* list_pp = (MMG5_pPoint*) malloc(parmesh->nprocs * sizeof(MMG5_pPoint));
  for (int k = 0; k < parmesh->nprocs; k++) list_pp[k] = (MMG5_pPoint) malloc(3*mesh->nt * sizeof(MMG5_Point));
  double** list_mm = (double**) malloc(parmesh->nprocs * sizeof(double*));
  for (int k = 0; k < parmesh->nprocs; k++) list_mm[k] = (double*) malloc(parmesh->met->size*3*mesh->nt * sizeof(double));
  int* list_deleted_triangles = (int*) malloc((mesh->nt+1) * sizeof(int));

  int ier, ieresult;
  int npt_global, npp_global;

  for (int i = 0; i < parmesh->nprocs; i++) {
    list_index_pt[i] = 0;
    list_index_pp[i] = 0;
  }
  for (int i = 0; i < mesh->nt; i++) list_deleted_triangles[i] = 0;

  // 1) Find the triangles to exchange
  int n = 0;
  for (int k = 1; k <= mesh->nt; k++) {
    pt = &mesh->tria[k];
    if ( !MG_EOK(pt) ) continue;

    if (partTab[k-1] == parmesh->myrank) continue;

    int send_proc = partTab[k-1];

    list_deleted_triangles[n++] = k;

    list_pt[send_proc][list_index_pt[send_proc]++] = mesh->tria[k];

    for (int i = 0; i <= 2; i++) {
      list_pp[send_proc][list_index_pp[send_proc]] = mesh->point[pt->v[i]];
      if (parmesh->met->size == 1) { // isotropic
        list_mm[send_proc][list_index_pp[send_proc]] = parmesh->met->m[pt->v[i]];
      }
      else { // anisotropic
        list_mm[send_proc][3*list_index_pp[send_proc]] = parmesh->met->m[3*pt->v[i]];
        list_mm[send_proc][3*list_index_pp[send_proc]+1] = parmesh->met->m[3*pt->v[i]+1];
        list_mm[send_proc][3*list_index_pp[send_proc]+2] = parmesh->met->m[3*pt->v[i]+2];
      }
      list_index_pp[send_proc]++;
    }
  }

  // 2) Exchange the triangles, points and metrics
  ier = PMMG2D_exchange( parmesh, mesh->nt, list_index_pt, list_pt, list_index_pp, list_pp, list_mm,
                         &pt_global, &pp_global, &mm_global, &npt_global, &npp_global );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error in the mpi exchanges during the ParMetis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  for (int k = 0; k < parmesh->nprocs; k++) {
    free(list_pt[k]);
    free(list_pp[k]);
    free(list_mm[k]);
  }
  free(list_pt);
  free(list_pp);
  free(list_mm);

  // 3) Remove the sent triangles from the local mesh
  int j = 0;
  for (int k = 1; k <= mesh->nt; k++ ) {
    if (list_deleted_triangles[j] == k) {
      j += 1;
      continue;
    }
    mesh->tria[k-j] = mesh->tria[k];
  }

  free(list_deleted_triangles);

  mesh->nt -= j;
  int nt_tmp = mesh->nt - j; // temporary number of triangles, used to check the contiguity (removing twice the triangles)

  PMMG2D_REALLOC(mesh, mesh->tria, mesh->nt+1, mesh->ntmax+1,
                 MMG5_Tria,"ParMetis decomposition realloc triangles",return 0);
  PMMG2D_REALLOC(mesh, mesh->point, mesh->np+1, mesh->npmax+1,
                 MMG5_Point,"ParMetis decomposition realloc points", return 0);
  PMMG2D_REALLOC(mesh, parmesh->met->m, parmesh->met->size*(mesh->np+1), parmesh->met->size*(mesh->npmax+1),
                 double, "ParMetis decomposition realloc metrics", return 0);

  // 4) Remove the points that do not belong to the mesh anymore
  ier = PMMG2D_remove_points( parmesh );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when removing the points during the ParMetis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 5) Add the triangles to their corresponding new process
  ier = PMMG2D_add_triangles( parmesh, pt_global, pp_global, mm_global, npt_global, npp_global);
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when adding the new triangles during the ParMetis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  PMMG2D_update( parmesh );

  // 6) Check the contiguity and correct the isolated triangles
  ier = PMMG2D_check_contiguity( parmesh, nt_tmp );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when checking the contiguity during the ParMetis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  PMMG2D_update( parmesh );

  // 7) Untag the interface parallel points and edges
  PMMG2D_untag_parallel( parmesh );

  // 8) Mark the interface nodes
  ier = PMMG2D_mark_parallel_interface_nodes( parmesh );
  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );
  if ( !ieresult ) {
    fprintf(stderr,"\n  ## Error when marking the interface nodes during the contiguity check in the ParMetis decomposition. Exit program.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // 9) Update the list of interface nodes
  PMMG2D_free_interface_nodes_list( parmesh );
  PMMG2D_fill_interface_nodes_list( parmesh );

  // Free the pointers
  PMMG2D_DEL_MEM(mesh, pt_global, MMG5_Tria, "Delete pt_global in the ParMetis decomposition");
  PMMG2D_DEL_MEM(mesh, pp_global, MMG5_Point, "Delete pp_global in the ParMetis decomposition");
  PMMG2D_DEL_MEM(mesh, mm_global, double, "Delete mm_global in the ParMetis decomposition");

  // 10) Rebuild the boundary edges
  if (mesh->edge) {
    PMMG2D_DEL_MEM(mesh, mesh->edge, MMG5_Edge, "Delete edges when decomposing the mesh using ParMetis");
    PMMG2D_REALLOC(mesh, mesh->edge, 1, 0, MMG5_Edge, "Realloc edge when decomposing the mesh using ParMetis (1)", return 0);
    mesh->na = 0;
  }

  for (int k = 1; k <= mesh->nt; k++ ) {
    pt = &mesh->tria[k];

    if ( !MG_EOK(pt) ) continue;

    for (int i = 0; i < 3; i++ ) {

      if (pt->tag[i] & MG_BDY) {
        mesh->na++;
        PMMG2D_REALLOC(mesh, mesh->edge, mesh->na+1, mesh->na,
                       MMG5_Edge, "Realloc edge when decomposing the mesh using ParMetis (2)", return 0);
        int i1 = MMG5_inxt2[i];
        int i2 = MMG5_inxt2[i1];
        mesh->edge[mesh->na].a = pt->v[i1];
        mesh->edge[mesh->na].b = pt->v[i2];
        mesh->edge[mesh->na].tag = pt->tag[i];
        mesh->edge[mesh->na].ref = mesh->na;
        pt->edg[i] = mesh->na;
      }
    }
  }

  return 1;
}


/**
 * \param parmesh pointer toward the parmesh structure.
 * \return 0 if it fails, 1 otherwise. 
 * 
 * Decompose the mesh in parallel using ParMetis to balance the loads.
 *
 */
int PMMG2D_ParMetis_decomposition( PMMG2D_pParMesh parmesh )
{
  MMG5_pMesh mesh;
  mesh = parmesh->mesh;
  int ier, ieresult;
  MPI_Status status;

  idx_t nparts = parmesh->nprocs;
  idx_t numflag = 0;
  idx_t ncon = 1;
  idx_t *vwgt = NULL;
  idx_t *adjwgt = NULL;
  idx_t wgtflag = 0;

  real_t *tpwgts = (real_t *)malloc(nparts * sizeof(real_t));
  real_t *ubvec = (real_t *)malloc(ncon * sizeof(real_t));
  idx_t options[METIS_NOPTIONS];

  for (int i = 0; i < nparts; i++) {
    tpwgts[i] = 1.0 / nparts; // Equal partition weights
  }
  for (int i = 0; i < ncon; i++) ubvec[i] = (1.+PMMG2D_LOAD_IMBALANCE); // Load imbalance tolerance

  METIS_SetDefaultOptions(options); // Default options
  options[1] = 0;

  idx_t edgecut; // Output: Number of edges cut by the partitioning
  idx_t *xadj = (idx_t *)malloc((mesh->nt+1) * sizeof(idx_t));
  idx_t *part = (idx_t *)malloc(mesh->nt * sizeof(idx_t)); // Output: Partition assignment

  // List of the triangles that are the neighbours
  int* nb_neighbours = malloc(mesh->nt*sizeof(int));
  int NB_TRIANGLES = 3; // Maximum number of neighbouring triangles associated to each triangle
  int **list_neighbours = malloc(mesh->nt*sizeof(int*));
  for (int i = 0; i < mesh->nt; i++) {
    nb_neighbours[i] = 0;
    list_neighbours[i] = malloc(NB_TRIANGLES*sizeof(int));
    for (int j = 0; j < NB_TRIANGLES; j++) list_neighbours[i][j] = 0;
  }

  if (!PMMG2D_find_neighbour_triangles(parmesh, nb_neighbours, list_neighbours )) {
    fprintf(stderr, "Error when finding the neighbouring triangles in ParMetis decomposition.\n");
    PMMG2D_CLEAN_AND_RETURN(parmesh,PMMG2D_STRONGFAILURE);
  }

  // Fill the points and edges structures
  int numEdges = 0;
  for (int k = 0; k < mesh->nt; k++) numEdges += nb_neighbours[k];
  idx_t *adjncy = (idx_t *)malloc(numEdges * sizeof(idx_t));

  xadj[0] = 0;
  for (int k = 0; k < mesh->nt; k++) {
    xadj[k+1] = xadj[k] + nb_neighbours[k];
    for (int i = 0; i < nb_neighbours[k]; i++) {
      adjncy[xadj[k]+i] = list_neighbours[k][i];
    }
  }

  int rcounts[parmesh->nprocs];
  MPI_Allgather(&mesh->nt, 1, MPI_INT, rcounts, 1, MPI_INT, parmesh->comm);
  idx_t vtxdist[nparts+1]; // Distribution of vertices among processes
  vtxdist[0] = 0;
  for (int k = 1; k <= nparts; k++) vtxdist[k] = vtxdist[k-1] + rcounts[k-1];

  // Call ParMETIS
#ifdef USE_PARMETIS
  ParMETIS_V3_PartKway(
    vtxdist, xadj, adjncy, vwgt, adjwgt,
    &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec,
    options, &edgecut, part, &parmesh->comm
  );
#endif

  for (int i = 0; i < mesh->nt; i++) free(list_neighbours[i]);
  free(list_neighbours);
  free(nb_neighbours);

  // Exchange the triangles to each new process
  ier = PMMG2D_ParMetis_exchange(parmesh,part);

  MPI_Allreduce( &ier, &ieresult, 1, MPI_INT, MPI_MIN, parmesh->comm );

  // Free memory
  free(tpwgts);
  free(ubvec);
  free(part);
  free(xadj);
  free(adjncy);

  return ieresult;

}
