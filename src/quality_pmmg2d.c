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

#include "parmmg2d.h"

typedef struct {
  double avlen,lmin,lmax;
  int    ned,amin,bmin,amax,bmax,nullEdge,hl[9];
  int    cpu_min,cpu_max;
} PMMG2D_lenStats;

/**
 * \param parmesh pointer to parmesh structure
 * \param opt PMMG_INQUA if called before the Mmg call, PMMG_OUTQUA otherwise
 * \param isCentral 1 for centralized mesh (no parallel communication), 0 for
 * distributed mesh
 *
 * \return 1 if success, 0 if fail;
 *
 * Print quality histogram among all group meshes and all processors
 */
int PMMG2D_qualhisto( PMMG2D_pParMesh parmesh, int isCentral )
{
  MMG5_pTria    pt;
  double        rap,rapmin,rapmax,rapavg,med,good;
  double        rapmin_result, rapmax_result, rapavg_result, med_result, good_result;
  int           i,ir,imax,his[PMMG2D_QUAL_HISSIZE],his_result[PMMG2D_QUAL_HISSIZE],rank_min,rank_send;
  int           k,iel,ok;
  int64_t       nex,nex_result,nt_result;
  static int8_t mmgWarn0 = 0;

  // Compute triangle quality
  for (k=1; k <= parmesh->mesh->nt; k++) {
    pt = &parmesh->mesh->tria[k];
    if( !MG_EOK(pt) )   continue;

    if ( !parmesh->met->m || parmesh->met->size == 1 ) {
      pt->qual = MMG2D_caltri_iso(parmesh->mesh, parmesh->met, pt);
    }
    else
      pt->qual = MMG2D_caltri_ani(parmesh->mesh, parmesh->met, pt);
  }
  if ( parmesh->info.imprim0 <= PMMG2D_VERB_VERSION ) return 1;

  rapmin  = 2.0;
  rapmax  = 0.0;
  rapavg  = med = good = 0.0;
  iel     = 0;

  for (k = 0; k < PMMG2D_QUAL_HISSIZE; k++)  his[k] = 0;

  nex = ok = 0;
  for (k = 1; k <= parmesh->mesh->nt; k++) {
    pt = &parmesh->mesh->tria[k];
    if( !MG_EOK(pt) ) {
      nex++;
      continue;
    }
    ok++;

    if ( (!mmgWarn0) && (MMG2D_quickcal(parmesh->mesh,pt) < 0.0) ) {
      mmgWarn0 = 1;
      fprintf(stderr,"  ## Warning: %s: at least 1 negative area\n",__func__);
    }

    if ( !parmesh->met->m || parmesh->met->size == 1 ) {
      rap = MMG2D_ALPHAD * MMG2D_caltri_iso(parmesh->mesh, parmesh->met, pt);
    }
    else
      rap = MMG2D_ALPHAD * MMG2D_caltri_ani(parmesh->mesh, parmesh->met, pt);

    if ( rap < rapmin ) {
      rapmin = rap;
      iel    = ok;
    }
    if ( rap > 0.5 )  med++;
    if ( rap > 0.12 ) good++;
    if ( rap < MMG2D_BADKAL )  parmesh->mesh->info.badkal = 1;
    rapavg += rap;
    rapmax  = MG_MAX(rapmax,rap);
    ir = MG_MIN(4,(int)(5.0*rap));
    his[ir] += 1;
  }

  if (isCentral) {
    nt_result = parmesh->mesh->nt;
    nex_result = nex;
    rapmin_result = rapmin;
    rapmax_result = rapmax;
    rapavg_result = rapavg;
    good_result = good;
    med_result = med;
    for ( i = 0; i < PMMG2D_QUAL_HISSIZE; ++i ) his_result[i] = his[i];
  }
  else {
    MPI_Reduce( &parmesh->mesh->nt, &nt_result, 1, MPI_INT64_T, MPI_SUM, 0, parmesh->comm );
    MPI_Reduce( &nex, &nex_result, 1, MPI_INT64_T, MPI_SUM, 0, parmesh->comm );
    MPI_Reduce( &rapmin, &rapmin_result, 1, MPI_DOUBLE, MPI_MIN, 0, parmesh->comm );
    MPI_Reduce( &rapmax, &rapmax_result, 1, MPI_DOUBLE, MPI_MAX, 0, parmesh->comm );
    MPI_Reduce( &rapavg, &rapavg_result, 1, MPI_DOUBLE, MPI_SUM, 0, parmesh->comm );
    MPI_Reduce( &good, &good_result, 1, MPI_DOUBLE, MPI_SUM, 0, parmesh->comm );
    MPI_Reduce( &med, &med_result, 1, MPI_DOUBLE, MPI_SUM, 0, parmesh->comm );
    MPI_Reduce( his, his_result, PMMG2D_QUAL_HISSIZE, MPI_INT, MPI_SUM, 0, parmesh->comm );

    // To find the real worst element
    rank_min = -1;
    if (abs(rapmin - rapmin_result) < 1.e-12) rank_min = parmesh->myrank;
    MPI_Allreduce( &rank_min, &rank_send, 1, MPI_INT, MPI_MAX, parmesh->comm );
    MPI_Bcast(&iel, 1, MPI_INT, rank_send, parmesh->comm);
  }

  if ( parmesh->myrank == 0 ) {

    fprintf(stdout,"\n  -- MESH QUALITY   %" PRId64 "\n",nt_result - nex_result);
    fprintf(stdout,"     BEST   %8.6f  AVRG.   %8.6f  WRST.   %8.6f (%d)\n",
          rapmax_result, rapavg_result / (nt_result - nex_result), rapmin_result, iel);
    if ( parmesh->nprocs > 1 && !isCentral) fprintf( stdout, "     PROC %d - \n", rank_send);

    // print histo
    fprintf(stdout,"     HISTOGRAMM:");
    fprintf(stdout,"  %6.2f %% > 0.12\n",100.0*(good_result/(float)(nt_result - nex_result)));
    if ( abs(parmesh->info.imprim) > 3 ) {
      fprintf(stdout,"                  %6.2f %% >  0.5\n",100.0*( med_result/(float)(nt_result - nex_result)));
      imax = MG_MIN(4,(int)(PMMG2D_QUAL_HISSIZE*rapmax_result));
      for (i = imax; i >= (int)(PMMG2D_QUAL_HISSIZE*rapmin_result); i--) {
        fprintf(stdout,"     %5.1f < Q < %5.1f   %7d   %6.2f %%\n",
                i/(1.*PMMG2D_QUAL_HISSIZE),i/(1.*PMMG2D_QUAL_HISSIZE)+0.2, his_result[i],100.*(his_result[i]/(float)(nt_result - nex_result)));
      }
    }
    if (!MMG5_minQualCheck(iel,rapmin_result,1.)) return 0;
  }

  return 1;

}

/**
 * \param parmesh pointer to parmesh structure
 *
 * \return 1 if success, 0 if fail;
 *
 * Resume edge length histo computed on each procs on the root processor
 *
 * \warning for now, only callable on centralized mesh
 *
 */
int PMMG2D_prilen( PMMG2D_pParMesh parmesh )
{
  MMG5_pTria pt;
  double dned, len;
  int l, ipa, ipb;
  static double bd[9] = {0.0, 0.3, 0.6, 0.7071, 0.9, 1.3, 1.4142, 2.0, 5.0};
  PMMG2D_lenStats lenStats;

  lenStats.avlen = 0.; // lavg
  lenStats.lmin = DBL_MAX;
  lenStats.lmax = 0.;
  lenStats.ned = 0; //navg
  lenStats.amin = lenStats.amax = lenStats.bmin = lenStats.bmax = 0;
  lenStats.nullEdge = 0;
  memset(lenStats.hl,0,9*sizeof(int));
  lenStats.cpu_min = lenStats.cpu_max = parmesh->myrank;

  if ( parmesh->met && parmesh->met->m ) {
  
    for (int k = 1; k <= parmesh->mesh->nt; k++) {
      pt = &parmesh->mesh->tria[k];
      if ( !MG_EOK(pt) )  continue;
  
      for (int ia = 0; ia < 3; ia++) {
        l = (&parmesh->mesh->adja[3*(k-1)+1])[ia];
        if ( l < 3*k )  continue;
  
        ipa = MMG2D_iare[ia][0];
        ipb = MMG2D_iare[ia][1];

        if ( !parmesh->met->m ) {
          len = MMG2D_lencurv_iso(parmesh->mesh,parmesh->met,pt->v[ipa],pt->v[ipb]);
        }
        else {
          len = MMG2D_lencurv_ani(parmesh->mesh,parmesh->met,pt->v[ipa],pt->v[ipb]);
        }
 
        lenStats.ned++;
        lenStats.avlen += len;
  
        // find largest, smallest edge
        if (len < lenStats.lmin) {
          lenStats.lmin  = len;
          lenStats.amin = pt->v[ipa];
          lenStats.bmin = pt->v[ipb];
        }
        if (len > lenStats.lmax) {
          lenStats.lmax  = len;
          lenStats.amax = pt->v[ipa];
          lenStats.bmax = pt->v[ipb];
        }
  
        // update histogram
        if (len < bd[3]) {
          if (len > bd[2])       lenStats.hl[2]++;
          else if (len > bd[1])  lenStats.hl[1]++;
          else                   lenStats.hl[0]++;
        }
        else if (len < bd[5]) {
          if (len > bd[4])       lenStats.hl[4]++;
          else if (len > bd[3])  lenStats.hl[3]++;
        }
        else if (len < bd[6])    lenStats.hl[5]++;
        else if (len < bd[7])    lenStats.hl[6]++;
        else if (len < bd[8])    lenStats.hl[7]++;
        else                     lenStats.hl[8]++;
      }
    }

  }

  if ( parmesh->myrank == parmesh->info.root ) {
    dned                  = (double)lenStats.ned;
    lenStats.avlen = lenStats.avlen / dned;

    fprintf(stdout,"\n  -- RESULTING EDGE LENGTHS (ROUGH EVAL.) %d \n",lenStats.ned);
    fprintf(stdout,"     AVERAGE LENGTH         %12.4f\n",lenStats.avlen);
    fprintf(stdout,"     SMALLEST EDGE LENGTH   %12.4f   %6d %6d",
            lenStats.lmin,lenStats.amin,lenStats.bmin);
    if ( parmesh->nprocs>1 ) {
      fprintf(stdout," (PROC %d)\n",lenStats.cpu_min);
    }
    else { fprintf(stdout,"\n"); }

    fprintf(stdout,"     LARGEST  EDGE LENGTH   %12.4f   %6d %6d",
            lenStats.lmax,lenStats.amax,lenStats.bmax);
    if ( parmesh->nprocs>1 ) {
      fprintf(stdout," (PROC %d)\n",lenStats.cpu_max);
    }
    else { fprintf(stdout,"\n"); }


    MMG5_displayLengthHisto_internal ( lenStats.ned,lenStats.amin,lenStats.bmin,
                                       lenStats.lmin,lenStats.amax,
                                       lenStats.bmax,lenStats.lmax,
                                       lenStats.nullEdge,bd,
                                       lenStats.hl,1,parmesh->info.imprim);
  }

  return 1;
}

