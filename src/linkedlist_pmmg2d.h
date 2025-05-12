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
 * \file linkedlist_pmmg2d.h
 * \brief linkedlist_pmmg2d.c header file
 * \author Fabien Salmon (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#ifndef LINKEDLIST_PMMG2D_H

#define LINKEDLIST_PMMG2D_H

#include "parmmg2d.h"

#define PMMG2D_LISTSIZE 30

/**
 * \struct PMMG2D_lnkdCell
 *
 * \brief Cell of a linked list. This cell allow to store 2 values and an index.
 *
 */
typedef struct {
  int     val1; /*!< first value of the cell */
  int     val2; /*!< second value of the cell */
  int     id;   /*!< position of the cell in the array */
  int     nxt;  /*!< position of the next cell in the linked list */
} PMMG2D_lnkdCell;

/**
 * \struct PMMG2D_lnkdVal
 *
 * \brief Cell of a linked list. This cell allow to store 1 value.
 *
 */
typedef struct {
  int     val; /*!< value of the cell */
  int     nxt;  /*!< position of the next cell in the linked list */
} PMMG2D_lnkdVal;

/**
 * \struct PMMG2D_lnkdList
 *
 * \brief Linked list of Cells of. Each cell can store 2 values.
 *
 */
typedef struct {
  int nitem; /*!< number of used cells in the list (=position to insert the next cell) */
  int nitem_max; /*!< maximal number of item in the list */
  int frst; /*!< position of the first cell of the linked list */
  int id;  /*!< ID of the linked list */
  PMMG2D_lnkdCell *item; /*!< array of cells */
} PMMG2D_cellLnkdList;

/**
 * \struct PMMG_lnkdList
 *
 * \brief Linked list of cells. Each cell can store 2 values.
 *
 */
typedef struct {
  int nitem; /*!< number of used cells in the list (=position to insert the next cell) */
  int nitem_max; /*!< maximal number of item in the list */
  int frst; /*!< position of the first cell of the linked list */
  int id;  /*!< ID of the linked list */
  PMMG2D_lnkdVal *item; /*!< array of Values */
} PMMG2D_valLnkdList;

int  PMMG2D_valLnkdListNew( PMMG2D_pParMesh parmesh, PMMG2D_valLnkdList *list,int,int );
int  PMMG2D_add_val2lnkdList( PMMG2D_pParMesh, PMMG2D_valLnkdList*,int);
int  PMMG2D_pop_val_lnkdList( PMMG2D_pParMesh, PMMG2D_valLnkdList*,int*);
int  PMMG2D_compare_valLnkdListLen (const void * a, const void * b);

#endif
