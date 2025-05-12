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
 * \file linkedlist_pmmg2d.c
 * \brief functions to manage a linked list
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "linkedlist_pmmg2d.h"

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_lnkdList structure.
 * \param nitem_max maximal number of cells of the linked list.
 *
 * \return 1 if success, 0 if fail.
 *
 * Initialisation of a linked list of cells.
 *
 */
int PMMG2D_valLnkdListNew( PMMG2D_pParMesh parmesh, PMMG2D_valLnkdList *list,int id,int nitem_max )
{

  /* adjust hash table params */
  list->nitem     = 0;
  list->nitem_max = nitem_max;
  list->id        = id;

  PMMG2D_MALLOC(parmesh,list->item,nitem_max,PMMG2D_lnkdVal,"linked list array",
                return 0);

  list->frst = PMMG2D_UNSET;

  return 1;
}

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG2D_lnkdList structure.
 * \param val1      value to add to the \a val1 field of the new cell.
 * \param val2      value to add to the \a val2 field of the new cell.
 *
 * \return 1 if we insert a new cell and 0 if fail.
 *
 * Non sorted insertion of a new cell with values \a val1 and \a val2
 * in the linked list \a list.
 *
 */
int PMMG2D_add_val2lnkdList( PMMG2D_pParMesh parmesh, PMMG2D_valLnkdList *list,
                             int val ) {
  PMMG2D_lnkdVal *cell;
  int             k,newsize;

  // Start from the first cell of the list

  if ( list->frst < 0 ) {
    // Add the first element
    list->frst    = 0;
    cell          = &list->item[list->nitem++];
    cell->val     = val;
    cell->nxt     = -1;

    return 1;
  }

  // Cell insertion at the end of the list
  if ( list->nitem >= list->nitem_max ) {
    newsize = MG_MAX((int)(1.2*list->nitem_max),list->nitem_max+1);
    PMMG2D_REALLOC(parmesh,list->item,newsize,
                   list->nitem_max,PMMG2D_lnkdVal,"linked list",return 0);
    list->nitem_max = newsize;
  }
  k             = list->nitem++;
  cell          = &list->item[k];
  cell->val     = val;
  cell->nxt     = -1;

  list->item[k-1].nxt = k;

  return 1;
}

/**
 * \param a  pointer toward a PMMG2D_lnkdList* structure.
 * \param b  pointer toward a PMMG2D_lnkdList* structure.
 *
 * \return 1 if a is greater than b, -1 if b is greater than 1, 0 if they are
 * equals.
 *
 * Compare lengths of 2 linked lists (can be used inside the qsort C fnuction)
 *
 */
int PMMG2D_compare_valLnkdListLen (const void * a, const void * b) {
  PMMG2D_valLnkdList *list1,*list2;

  list1 = *(PMMG2D_valLnkdList**)a;
  list2 = *(PMMG2D_valLnkdList**)b;

  if ( list1->nitem > list2->nitem ) return 1;

  if ( list1->nitem < list2->nitem ) return -1;

  return 0;
}

/**
 * \param parmesh   pointer toward the parmesh structure.
 * \param list      pointer toward a PMMG_valLnkdList structure.
 * \param val       first list value (removed)
 *
 * \return 1
 *
 * Remove the first item of a list and fill \a val with this removed value.
 *
 */
int PMMG2D_pop_val_lnkdList( PMMG2D_pParMesh parmesh, PMMG2D_valLnkdList *list,
                             int *val ) {
  PMMG2D_lnkdVal *cell;

  // Get first cell
  cell = &list->item[list->frst];
  *val = cell->val;

  // Pop cell from the head of the list
  list->frst = cell->nxt;
  list->nitem--;

  return 1;
}

