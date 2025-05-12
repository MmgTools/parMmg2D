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
 * \file inout_pmmg2d.c
 * \brief io for the parmmg2d software
 * \author Fabien Salmon (Inria)
 * \version 1
 * \date 04 2024
 * \copyright GNU Lesser General Public License.
 *
 * input/outputs for parmmg2d.
 *
 */

#include "parmmg2d.h"

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the mesh from.
 * \param ftmin format of the input file.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load a centralized mesh. 
 *
 */
int PMMG2D_loadMesh_centralized(PMMG2D_pParMesh parmesh, const char *filename, int fmtin) {
  int        ier;
  const char *data;

  if ( parmesh->myrank != parmesh->info.root ) return 1;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( parmesh->mesh->info.imprim == parmesh->info.mmg_imprim );
  parmesh->mesh->info.imprim = MG_MAX ( parmesh->info.imprim, parmesh->mesh->info.imprim );

  if ( filename ) {
    data = filename;
  }
  else if ( parmesh->meshin ) {
    data = parmesh->meshin;
  }
  else if ( parmesh->mesh->namein ) {
    data = parmesh->mesh->namein;
  }
  else {
    data = NULL;
  }

  switch ( fmtin ) {
  case ( MMG5_FMT_GmshASCII ): case ( MMG5_FMT_GmshBinary ):
    ier = MMG2D_loadMshMesh(parmesh->mesh,parmesh->met,data);
    break;

  case ( MMG5_FMT_VtkVtp ):
    ier = MMG2D_loadVtpMesh(parmesh->mesh,parmesh->met,data);
    break;

  case ( MMG5_FMT_VtkVtu ):
    ier = MMG2D_loadVtuMesh(parmesh->mesh,parmesh->met,data);
    break;

  case ( MMG5_FMT_VtkVtk ):
    ier = MMG2D_loadVtkMesh(parmesh->mesh,parmesh->met,data);
    break;

  case ( MMG5_FMT_MeditASCII ): case ( MMG5_FMT_MeditBinary ):
    ier = MMG2D_loadMesh(parmesh->mesh,data);
    break;

  default:
    if ( parmesh->myrank == parmesh->info.root ) {
      fprintf(stderr,"  ** I/O AT FORMAT %s NOT IMPLEMENTED.\n",MMG5_Get_formatName(fmtin) );
    }
    ier = 0;
  }

  // Restore the mmg verbosity to its initial value
  parmesh->mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the metric from.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load a centralized metric. 
 *
 */
int PMMG2D_loadMet_centralized(PMMG2D_pParMesh parmesh, const char *filename) {
  int        ier;
  const char *data;

  if ( parmesh->myrank != parmesh->info.root ) return 1;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( parmesh->mesh->info.imprim == parmesh->info.mmg_imprim );
  parmesh->mesh->info.imprim = MG_MAX ( parmesh->info.imprim, parmesh->mesh->info.imprim );

  if ( filename ) {
    data = filename;
  }
  else if ( parmesh->metin ) {
    data = parmesh->metin;
  }
  else if ( parmesh->met->namein ) {
    data = parmesh->met->namein;
  }
  else {
    data = NULL;
  }
  ier = MMG2D_loadSol(parmesh->mesh, parmesh->met, data);

  /* Restore the mmg verbosity to its initial value */
  parmesh->mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the mesh from.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save a centralized mesh. 
 *
 */
int PMMG2D_saveMesh_centralized(PMMG2D_pParMesh parmesh, const char *filename) {
  int        ier;

  if ( parmesh->myrank!=parmesh->info.root ) return 1;

  // Set mmg verbosity to the max between the Parmmg2d verbosity and the mmg verbosity
  assert ( parmesh->mesh->info.imprim == parmesh->info.mmg_imprim );
  parmesh->mesh->info.imprim = MG_MAX ( parmesh->info.imprim, parmesh->mesh->info.imprim );

  if ( filename && *filename ) {
    ier = MMG2D_saveMesh(parmesh->mesh,filename);
  }
  else {
    ier = MMG2D_saveMesh(parmesh->mesh,parmesh->meshout);
  }

  /* Restore the mmg verbosity to its initial value */
  parmesh->mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the metric from.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save a centralized metric. 
 *
 */
int PMMG2D_saveMet_centralized(PMMG2D_pParMesh parmesh, const char *filename) {
  int        ier;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( parmesh->mesh->info.imprim == parmesh->info.mmg_imprim );
  parmesh->mesh->info.imprim = MG_MAX ( parmesh->info.imprim, parmesh->mesh->info.imprim );

  if ( filename && *filename ) {
    ier =  MMG2D_saveSol(parmesh->mesh,parmesh->met,filename);
  }
  else {
    ier =  MMG2D_saveSol(parmesh->mesh,parmesh->met,parmesh->metout);
  }

  /* Restore the mmg verbosity to its initial value */
  parmesh->mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save a centralized mesh with the metrics in a vtk file. 
 *
 */
int PMMG2D_saveMeshandMet_centralized_vtk(PMMG2D_pParMesh parmesh) {
  int k;
  MMG5_pTria pt;
  MMG5_pEdge pe;

  if ( parmesh->myrank!=parmesh->info.root ) {
    return 1;
  }

  FILE *fp = fopen(parmesh->meshout, "w");
  fprintf(fp, "# vtk DataFile Version 3.0\nMesh with metrics at points\nASCII\nDATASET UNSTRUCTURED_GRID\n");
  fprintf(fp, "POINTS %d double\n", parmesh->mesh->np);

  // Vertices
  for (k = 1; k < parmesh->mesh->np+1; k++) {
    fprintf(fp, "%f %f 0.\n", parmesh->mesh->point[k].c[0], parmesh->mesh->point[k].c[1]);
  }

  fprintf(fp, "\nCELLS %d %d\n", parmesh->mesh->nt + parmesh->mesh->na, 4*parmesh->mesh->nt + 3*parmesh->mesh->na);

  // Triangles
  for (k = 1; k < parmesh->mesh->nt+1; k++) {
    pt = &parmesh->mesh->tria[k];
    fprintf(fp, "3 %d %d %d \n",pt->v[0]-1, pt->v[1]-1, pt->v[2]-1);
  }

  // Edges
  for (k = 1; k < parmesh->mesh->na+1; k++) {
    pe = &parmesh->mesh->edge[k];
    fprintf(fp, "2 %d %d \n",pe->a-1, pe->b-1);
  }

  // Cell types
  fprintf(fp, "\nCELL_TYPES %d\n", parmesh->mesh->nt + parmesh->mesh->na);
  for (int k = 1; k < parmesh->mesh->nt+1; k++) fprintf(fp, "5\n");
  for (int k = 1; k < parmesh->mesh->na+1; k++) fprintf(fp, "3\n");

  if (!parmesh->met && !parmesh->met->m) {
    fclose(fp);
    return 1;
  }

  // Metrics
  fprintf(fp, "\nPOINT_DATA %d\n", parmesh->mesh->np);
  fprintf(fp, "VECTORS Metric double \n");
  if (parmesh->met->size == 1) {
    for (int n = 1; n < parmesh->mesh->np+1; n++) {
      fprintf(fp, "%f %f %f\n", parmesh->met->m[n], 0., 0.);
    }
  }
  else {
    for (int n = 1; n < parmesh->mesh->np+1; n++) {
      fprintf(fp, "%f %f %f\n", parmesh->met->m[3*n], parmesh->met->m[3*n+1], parmesh->met->m[3*n+2]);
    }
  }

  fclose(fp);

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save a distributed mesh with the metrics in several mesh files.
 *
 */
int PMMG2D_saveMeshandMet_distributed(PMMG2D_pParMesh parmesh) {
  int ier = MMG2D_saveMesh(parmesh->mesh,parmesh->meshout);
  ier =  MMG2D_saveSol(parmesh->mesh,parmesh->met,parmesh->metout);
  return ier;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Save a distributed mesh with the metrics in several vtk files.
 *
 */
int PMMG2D_saveMeshandMet_distributed_vtk(PMMG2D_pParMesh parmesh) {
  int k;
  MMG5_pTria pt;
  MMG5_pEdge pe;

  FILE *fp = fopen(parmesh->meshout, "w");
  fprintf(fp, "# vtk DataFile Version 3.0\nMesh with metrics at points\nASCII\nDATASET UNSTRUCTURED_GRID\n");
  fprintf(fp, "POINTS %d double\n", parmesh->mesh->np);

  // Vertices
  for (k = 1; k < parmesh->mesh->np+1; k++) {
    fprintf(fp, "%f %f 0.\n", parmesh->mesh->point[k].c[0], parmesh->mesh->point[k].c[1]);
  }

  fprintf(fp, "\nCELLS %d %d\n", parmesh->mesh->nt + parmesh->mesh->na, 4*parmesh->mesh->nt + 3*parmesh->mesh->na);

  // Triangles
  for (k = 1; k < parmesh->mesh->nt+1; k++) {
    pt = &parmesh->mesh->tria[k];
    fprintf(fp, "3 %d %d %d \n",pt->v[0]-1, pt->v[1]-1, pt->v[2]-1);
  }

  // Edges
  for (k = 1; k < parmesh->mesh->na+1; k++) {
    pe = &parmesh->mesh->edge[k];
    fprintf(fp, "2 %d %d \n",pe->a-1, pe->b-1);
  }

  // Cell types
  fprintf(fp, "\nCELL_TYPES %d\n", parmesh->mesh->nt + parmesh->mesh->na);
  for (int k = 1; k < parmesh->mesh->nt+1; k++) fprintf(fp, "5\n");
  for (int k = 1; k < parmesh->mesh->na+1; k++) fprintf(fp, "3\n");

  if (!parmesh->met && !parmesh->met->m) {
    fclose(fp);
    return 1;
  }

  // Metrics
  fprintf(fp, "\nPOINT_DATA %d\n", parmesh->mesh->np);
  fprintf(fp, "VECTORS Metric double \n");
  if (parmesh->met->size == 1) {
    for (int n = 1; n < parmesh->mesh->np+1; n++) {
      fprintf(fp, "%f %f %f\n", parmesh->met->m[n], 0., 0.);
    }
  }
  else {
    for (int n = 1; n < parmesh->mesh->np+1; n++) {
      fprintf(fp, "%f %f %f\n", parmesh->met->m[3*n], parmesh->met->m[3*n+1], parmesh->met->m[3*n+2]);
    }
  }

  fclose(fp);

  return 1;
}

/**
 * \param n integer for which we want to know the number of digits
 *
 * \return the number of digits of n.
 *
 */
static inline
int PMMG2D_count_digits(int n) {

  int count = 0;
  while (n != 0) {
    n /= 10;
    ++count;
  }

  return count;
}

/**
 * \param parmesh pointer toward the parmesh structure
 * \param endame string to allocate to store the final filename
 * \param initname initial file name in which we want to insert rank index
 * \param ASCIIext extension to search for ASCII format
 * \param binext extension to search for binary format
 *
 * Allocate the endname string and copy the initname string with the mpir rank
 * index before the file extension.
 *
 */
static inline
void PMMG2D_insert_rankIndex(PMMG2D_pParMesh parmesh,char **endname,const char *initname,
                             char *ASCIIext, char *binext) {
  int    lenmax;
  int8_t fmt;
  char   *ptr;

  lenmax = PMMG2D_count_digits ( parmesh->nprocs );

  /* Check for pointer validity */
  if ( (!endname) || (!initname) ) {
    return;
  }

  MMG5_SAFE_CALLOC(*endname,strlen(initname)+lenmax+7,char,return);

  strcpy(*endname,initname);

  ptr = strstr(*endname,".meshb");

  fmt = 0; /* noext */
  if( ptr ) {
    *ptr = '\0';
    fmt = 1; /* binary */
  }
  else {
    ptr = strstr(*endname,".mesh");
    if( ptr ) {
      *ptr = '\0';
      fmt = 2; /* ASCII */
    }
  }
  sprintf(*endname, "%s.%d", *endname, parmesh->myrank );

  if ( fmt == 1 ) {
    strcat ( *endname, binext );
  }
  else if ( fmt == 2 ) {
    strcat ( *endname, ASCIIext );
  }

  return;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the mesh from.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load a distributed mesh in Medit format. 
 * The rank index is inserted in the input file name.
 *
 */
int PMMG2D_loadMesh_distributed(PMMG2D_pParMesh parmesh, const char *filename) {
  MMG5_pMesh  mesh;
  int         ier;
  char*       data = NULL;

  mesh = parmesh->mesh;

  // Add rank index to mesh name
  if ( filename ) {
    PMMG2D_insert_rankIndex(parmesh, &data, filename, ".mesh", ".meshb");
  }
  else if ( parmesh->meshin ) {
    PMMG2D_insert_rankIndex(parmesh, &data, parmesh->meshin, ".mesh", ".meshb");
  }
  else if ( mesh->namein ) {
    PMMG2D_insert_rankIndex(parmesh, &data, mesh->namein, ".mesh", ".meshb");
  }

  // Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier = MMG2D_loadMesh(mesh,data);

  // Restore the mmg verbosity to its initial value
  mesh->info.imprim = parmesh->info.mmg_imprim;

  if ( ier < 1 ) {
    MMG5_SAFE_FREE(data);
    return ier;
  }

  MMG5_SAFE_FREE(data);

  if ( 1 != ier ) return 0;

  return 1;
}

/**
 * \param parmesh pointer toward the parmesh structure.
 * \param filename name of the file to load the metric from.
 *
 * \return 0 if fail, 1 otherwise
 *
 * Load a distributed metric in Medit format. 
 * The rank index is inserted in the input file name.
 *
 */
int PMMG2D_loadMet_distributed(PMMG2D_pParMesh parmesh, const char *filename) {
  MMG5_pMesh mesh;
  MMG5_pSol  met;
  int        ier;
  char       *data = NULL;

  mesh = parmesh->mesh;
  met  = parmesh->met;

  // Add rank index to mesh name
  if ( filename ) {
    PMMG2D_insert_rankIndex(parmesh, &data, filename, ".sol", ".sol");
  }
  else if ( parmesh->metin ) {
    PMMG2D_insert_rankIndex(parmesh, &data, parmesh->metin, ".sol", ".sol");
  }
  else if ( met->namein ) {
    PMMG2D_insert_rankIndex(parmesh, &data, met->namein, ".sol", ".sol");
  }
  else if ( parmesh->meshin ) {
    PMMG2D_insert_rankIndex(parmesh, &data, parmesh->meshin, ".mesh", ".meshb");
  }
  else if ( mesh->namein ) {
    PMMG2D_insert_rankIndex(parmesh, &data, mesh->namein, ".mesh", ".meshb");
  }

  // Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity
  assert ( mesh->info.imprim == parmesh->info.mmg_imprim );
  mesh->info.imprim = MG_MAX ( parmesh->info.imprim, mesh->info.imprim );

  ier = MMG2D_loadSol(mesh,met,data);

  // Restore the mmg verbosity to its initial value 
  mesh->info.imprim = parmesh->info.mmg_imprim;

  MMG5_SAFE_FREE(data);

  return ier;
}
