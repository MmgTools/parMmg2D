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
 * \file libparmmg2d_tools.c
 * \brief C API functions definitions for PARMMG2D library.
 * \copyright GNU Lesser General Public License.
 *
 * C API for PARMMG2D library.
 *
 */
#include "parmmg2d.h"

/*! Helper macro used only in this file: copies the contents of fromV[fromC]
 *  to toV[toC] updates toC */
#define ARGV_APPEND(parmesh,fromV,toV,fromC,toC,msg,on_failure)   do {  \
    PMMG2D_MALLOC(parmesh, toV[ toC ], strlen( fromV[ fromC ] ) + 1, char, msg, \
                on_failure);                                            \
    memcpy( toV[ toC ], fromV[ fromC ], (strlen( fromV[ fromC ] ) + 1)*sizeof(char) ); \
    ++toC;                                                              \
  }while(0)

/**
 * \param parmesh pointer to pmmg structure
 * \param mmgArgv pointer to argv like buffer
 * \param mmgArgc pointer to argc like buffer
 * \param argc    actual argc value
 *
 * Free the allocations of the custom created argc/argv wrapper that is passed
 * to mmg to parse the command line options
 */
static void
PMMG2D_argv_cleanup( PMMG2D_pParMesh parmesh, char **mmgArgv, int mmgArgc, int argc )
{
  for (int i = 0; i < mmgArgc; ++i )
    PMMG2D_DEL_MEM(parmesh, mmgArgv[i],char, "Deallocating mmgargv[i]: " );
  PMMG2D_DEL_MEM(parmesh, mmgArgv,char*, "Deallocating mmgargv: " );
}

int PMMG2D_parsar( int argc, char *argv[], PMMG2D_pParMesh parmesh )
{
  int        val,i  = 0;
  int        ret_val = 1;
  int        mmgArgc = 0;
  char**     mmgArgv = NULL;

  /** Parse arguments specific to parMmg2d then add to mmgArgv the mmg arguments
   * and call the mmg2d parser. */
  for ( i = 1; i < argc; ++i ) {
    if ( !strcmp( argv[ i ],"-val" ) ) {
      RUN_ON_ROOT_AND_BCAST( PMMG2D_defaultValues(parmesh),0,
                             parmesh->myrank,ret_val=0; goto fail_mmgargv);
      ret_val = 0;
      goto fail_mmgargv;
    }
    else if ( ( !strcmp( argv[ i ],"-?" ) ) || ( !strcmp( argv[ i ],"-h" ) ) ) {
      RUN_ON_ROOT_AND_BCAST( PMMG2D_usage(parmesh, argv[0]),0,
                             parmesh->myrank,ret_val=0; goto fail_mmgargv);
      ret_val = 0;
      goto fail_mmgargv;
    }
  }

  /* Create a new set of argc/argv variables adding only the cl options that
     mmg has to process
     Overallocating as they are at most argc. Trying to avoid the overallocation
     is not worth any effort, these are ~kb */
  PMMG2D_MALLOC(parmesh, mmgArgv, argc, char*, " copy of argv for mmg: ",
                ret_val = 0; goto fail_mmgargv);

  /* First argument is always argv[0] ie prog name */
  i = 0;
  ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc, " mmgArgv[0] for mmg: ",
              ret_val = 0; goto fail_proc);

  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch( argv[i][1] ) {
      case 'c':
        if ( !strcmp(argv[i],"-centralized-output") ) {
          /* force centralized output: only relevant using medit distributed
           * input or library call */
          if ( !PMMG2D_Set_iparameter(parmesh,PMMG2D_IPARAM_distributedOutput,0) )  {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;
      // field not implemented yet
      /*case 'f':
        if ( !strcmp(argv[i],"-field") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-' ) {
            if ( ! PMMG2D_Set_inputSolsName(parmesh,argv[i]) ) {
              RUN_ON_ROOT_AND_BCAST( PMMG2D_usage(parmesh, argv[0]),0,
                                     parmesh->myrank,ret_val=0; goto fail_mmgargv);
              ret_val = 0;
              goto fail_mmgargv;
            }
          }
          else {
            RUN_ON_ROOT_AND_BCAST( PMMG2D_usage(parmesh, argv[0]),0,
                                   parmesh->myrank,ret_val=0; goto fail_mmgargv);
            ret_val = 0;
            goto fail_mmgargv;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;*/
      case 'h':
        if ( !strcmp(argv[i],"-hmin") && ++i < argc ) {
          if ( !PMMG2D_Set_dparameter(parmesh,PMMG2D_DPARAM_hmin,atof(argv[i])) ) {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-hmax") && ++i < argc ) {
          if ( !PMMG2D_Set_dparameter(parmesh,PMMG2D_DPARAM_hmax,atof(argv[i])) ) {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'i':
        if ( 0 == strncmp( argv[i], "-isotropic", 9 ) ) {
          parmesh->mesh->info.isotropic = 1;
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'l':
        if ( !strcmp(argv[i],"-limit-angle") && ++i < argc ) {
          if (!MMG2D_Set_dparameter(parmesh->mesh, parmesh->met, MMG2D_DPARAM_limit_angle, atof(argv[i])) ) {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'm':
        if ( !strcmp(argv[i],"-mmg-v") ) {

          /* Mmg verbosity */
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ||
                 (argv[i][0]=='-' && isdigit(argv[i][1])) ) {
              val = atoi(argv[i]);

              if ( !PMMG2D_Set_iparameter(parmesh,PMMG2D_IPARAM_mmgVerbose,val) ) {
                ret_val = 0;
                goto fail_proc;
              }
            }
            else {
              i--;
            }
          }
          else {
            fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
            ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-mmg-d") ) {
          if ( !PMMG2D_Set_iparameter(parmesh,PMMG2D_IPARAM_mmgDebug,val) ) {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-m") ) {
          /* memory */
          if ( ++i < argc && isdigit( argv[i][0] ) ) {
            if ( ( atoi(argv[ i ]) > MMG5_memSize() ) || ( atoi(argv[ i ]) < 0 ) ) {
              fprintf( stderr,
                       "\nErroneous mem size requested (%s)\n",argv[i] );
              ret_val = 0;
              goto fail_proc;
            }
            else {
              parmesh->info.mem = atoi( argv[i] );
              PMMG2D_parmesh_SetMemGloMax( parmesh );
            }
            PMMG2D_parmesh_SetMemMax( parmesh );
          } else {
            fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          /* else : what happens with -met option... to treat */
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'n':  /* number of adaptation iterations */
        if ( ( 0 == strncmp( argv[i], "-niter", 5 ) ) && ( ( i + 1 ) < argc ) ) {
          ++i;
          if ( isdigit( argv[i][0] ) && ( atoi( argv[i] ) >= 0 ) ) {
            parmesh->niter = atoi( argv[i] );
          } else {
            parmesh->niter = PMMG2D_NITER;
            fprintf( stderr, "\nWrong number of adaptation iterations (%s).\n",argv[i]);

            ret_val = 0;
            goto fail_proc;
          }
        } else if ( ( 0 == strncmp( argv[i], "-nlayers", 5 ) ) && ( ( i + 1 ) < argc ) ) {
          ++i;
          if ( isdigit( argv[i][0] ) && ( atoi( argv[i] ) > 0 ) ) {
            parmesh->info.ifc_layers = atoi( argv[i] );
          } else {
            parmesh->info.ifc_layers = PMMG2D_MVIFCS_NLAYERS;
            fprintf( stderr, "\nWrong number of layers for interface displacement (%s).\n",argv[i]);

            ret_val = 0;
            goto fail_proc;
          }
        } else if ( ( 0 == strncmp( argv[i], "-n_interp_layers", 12 ) ) && ( ( i + 1 ) < argc ) ) {
          ++i;
          if ( isdigit( argv[i][0] ) && ( atoi( argv[i] ) > 0 ) ) {
            parmesh->info.interp_layers = atoi( argv[i] );
          } else {
            parmesh->info.interp_layers = PMMG2D_INTERP_NLAYERS;
            fprintf( stderr, "\nWrong number of layers for interpolation (%s).\n",argv[i]);

            ret_val = 0;
            goto fail_proc;
          }
        } else if ( 0 == strncmp( argv[i], "-nobalance", 9 ) ) {
          parmesh->info.nobalancing = MMG5_ON;
        } else if ( 0 == strncmp( argv[i], "-nofem", 5 ) ) {
          if ( !PMMG2D_Set_iparameter(parmesh,PMMG2D_IPARAM_nofem,1) )  {
            ret_val = 0;
            goto fail_proc;
          }
        } else if ( 0 == strncmp( argv[i], "-noout", 5 ) ) {
          parmesh->info.fmtout = PMMG2D_UNSET;
        } else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'o':
        if ( 0 == strncmp( argv[i], "-optim-interp", 12 ) ) {
          parmesh->info.optim_interp = 1;
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'd':
        if ( !strcmp(argv[i],"-distributed-output") ) {
          /* force distributed output: only relevant using medit centralized
           * input or library call */
          if ( !PMMG2D_Set_iparameter(parmesh,PMMG2D_IPARAM_distributedOutput,1) )  {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else if ( !strcmp(argv[i],"-d") ) {
          /* debug */
          if ( !PMMG2D_Set_iparameter(parmesh,PMMG2D_IPARAM_debug,1) )  {
            ret_val = 0;
            goto fail_proc;
          }
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 'r':
        if ( !strcmp(argv[i],"-ratio-load-balance") ) {

          /* Load balancing tolerance for output redecomposition */
          if ( ++i < argc ) {
            if ( isdigit(argv[i][0]) ) {
              val = atof(argv[i]);

              if ( !PMMG2D_Set_dparameter(parmesh,PMMG2D_DPARAM_load_balance,val) ) {
                ret_val = 0;
                goto fail_proc;
              }
            }
            else {
              fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
              ret_val = 0;
              goto fail_proc;
            }
          }
          else {
            fprintf( stderr, "\nMissing argument option %c\n", argv[i-1][1] );
            ret_val = 0;
            goto fail_proc;
          }
        }

        if ( !strcmp(argv[i],"-rn") ) {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;

      case 's':
        if ( 0 == strncmp( argv[i], "-surf", 4 ) ) {
          parmesh->mesh->info.nosurf = 0;
        }
        else {
          ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                      " adding to mmgArgv for mmg: ",
                      ret_val = 0; goto fail_proc );
        }
        break;
      case 'v':  /* verbosity */
        if ( ++i < argc ) {
          if ( isdigit(argv[i][0]) ||
               (argv[i][0]=='-' && isdigit(argv[i][1])) ) {
            if ( !PMMG2D_Set_iparameter(parmesh,PMMG2D_IPARAM_verbose,atoi(argv[i])) ) {
              ret_val = 0;
              goto fail_proc;
            }
          }
          else
            i--;
        }
        else {
          fprintf(stderr,"\nMissing argument option %c\n",argv[i-1][1]);
          ret_val = 0;
          goto fail_proc;
        }
        break;

      default:
        ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                    " adding to mmgArgv for mmg: ",
                    ret_val = 0; goto fail_proc);

        break;
      }
    } else {
      ARGV_APPEND(parmesh, argv, mmgArgv, i, mmgArgc,
                  " adding to mmgArgv for mmg: ",
                  ret_val = 0; goto fail_proc);
    }
    ++i;
  }

  // parmmg finished parsing arguments, the rest will be handled by mmg2d
  if ( MMG2D_parsar( mmgArgc, mmgArgv, parmesh->mesh, parmesh->met, NULL ) != 1 ) {
    ret_val = 0;
    goto fail_proc;
  }

  if( parmesh->mesh->info.opnbdy ) {
    fprintf(stderr," ## Warning: Surface adaptation not supported with opnbdy."
        "\nSetting nosurf on.\n");
    if ( !MMG2D_Set_iparameter(parmesh->mesh,NULL,MMG2D_IPARAM_nosurf,1) ) return 0;
  }

  /* Store mesh names into the parmesh if needed */
  if ( !parmesh->meshin ) {
    assert ( parmesh->mesh->namein );
    PMMG2D_Set_name(parmesh,&parmesh->meshin,
                    parmesh->mesh->namein,"mesh.mesh");
  }
  if ( !parmesh->meshout ) {
    assert ( parmesh->mesh->nameout );
    PMMG2D_Set_name(parmesh,&parmesh->meshout,
                    parmesh->mesh->nameout,"mesh.o.mesh");
  }
  if ( (!parmesh->metin) && parmesh->met && parmesh->met->namein ) {
    PMMG2D_Set_name(parmesh,&parmesh->metin,
                    parmesh->met->namein,"mesh.sol");
  }
  if ( (!parmesh->metout) && parmesh->met && parmesh->met->nameout ) {
    PMMG2D_Set_name(parmesh,&parmesh->metout,
                    parmesh->met->nameout,"mesh.o.sol");
  }
  // level set and displacement not supported yet
  if ( (!parmesh->lsin) && parmesh->ls && parmesh->ls->namein ) {
    PMMG2D_Set_name(parmesh,&parmesh->lsin,
                    parmesh->ls->namein,"mesh.sol");
  }
  if ( (!parmesh->dispin) && parmesh->disp && parmesh->disp->namein ) {
    PMMG2D_Set_name(parmesh,&parmesh->dispin,
                    parmesh->disp->namein,"mesh.sol");
  }

fail_proc:
  PMMG2D_argv_cleanup( parmesh, mmgArgv, mmgArgc, argc );
fail_mmgargv:
  return ret_val;
}

#undef ARGV_APPEND

int PMMG2D_usage( PMMG2D_pParMesh parmesh, char * const prog )
{
  if ( !parmesh->myrank ) {
    fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",
            prog);

    fprintf(stdout,"\n** Generic options :\n");
    fprintf(stdout,"-h         Print this message\n");
    fprintf(stdout,"-v [n]     Tune ParMmg2d level of verbosity, [-10..10]\n");
    fprintf(stdout,"-mmg-v [n] Tune Mmg level of verbosity, [-10..10]\n");
    fprintf(stdout,"-m [n]     Set maximal memory size to n Mbytes\n");
    fprintf(stdout,"-d         Turn on debug mode for ParMmg2d\n");
    fprintf(stdout,"-mmg-d     Turn on debug mode for Mmg\n");
    fprintf(stdout,"-val       Print the default parameters values\n");

    fprintf(stdout,"\n**  File specifications\n");
    fprintf(stdout,"-in    file  input triangulation\n");
    fprintf(stdout,"-out   file  output triangulation\n");
    fprintf(stdout,"-sol   file  load metric file\n");
    fprintf(stdout,"-noout       do not write output triangulation\n");

    fprintf(stdout,"\n**  Parameters\n");
    fprintf(stdout,"-niter               val  number of remeshing iterations\n");
    fprintf(stdout,"-nlayers             val  number of layers for interface displacement\n");
    fprintf(stdout,"-n_interp_layers     val  number of layers for interpolation\n");
    fprintf(stdout,"-nobalance                switch off load balancing of the output mesh\n");
    fprintf(stdout,"-ratio-load-balance  val  load balancing tolerance for output redecomposition\n");
    fprintf(stdout,"-optim-interp             switch on the optimisation of interpolation between meshes\n");

    fprintf(stdout,"-hmin         val  minimal mesh size\n");
    fprintf(stdout,"-hmax         val  maximal mesh size\n");
    fprintf(stdout,"-hsiz         val  constant mesh size\n");
    fprintf(stdout,"-hgrad        val  control gradation\n");
    fprintf(stdout,"-hgradreq     val  control gradation from required entities\n");
    fprintf(stdout,"-A                 enable anisotropy (without metric file).\n");

#ifdef USE_SCOTCH
    fprintf(stdout,"-rn [n]            Turn on or off the renumbering using SCOTCH [1/0] \n");
#endif
    fprintf(stdout,"\n");

    fprintf(stdout,"-nofem       do not force Mmg to create a finite element mesh \n");
    fprintf(stdout,"-optim       mesh optimization\n");
    fprintf(stdout,"-noinsert    no point insertion/deletion \n");
    fprintf(stdout,"-noswap      no edge or face flipping\n");
    fprintf(stdout,"-nomove      no point relocation\n");
    fprintf(stdout,"-nosurf      no surface modifications\n");
    fprintf(stdout,"\n\n");

  }

  return 1;
}

int PMMG2D_defaultValues( PMMG2D_pParMesh parmesh )
{
  int ier = 1;

  if ( !parmesh->myrank ) {
    fprintf(stdout,"\n\n");
    fprintf(stdout,"  --- ParMMG2D ---\n");
    fprintf(stdout,"default parameter values:\n");
    fprintf(stdout,"\n** Generic options\n");
    fprintf(stdout,"verbosity                 (-v)      : %d\n",
            parmesh->info.imprim);

    fprintf(stdout,"maximal memory size       (-m)      : %zu MB\n",
            parmesh->memGloMax/MMG5_MILLION);
    fprintf(stdout,"\n** Parameters\n");
    fprintf( stdout,"# of remeshing iterations (-niter)        : %d\n",parmesh->niter);
    fprintf( stdout,"# of layers for interface displacement (-nlayers) : %d\n",PMMG2D_MVIFCS_NLAYERS);

#ifdef USE_SCOTCH
    fprintf(stdout,"SCOTCH renumbering                  : enabled\n");
#else
    fprintf(stdout,"SCOTCH renumbering                  : disabled\n");
#endif

    if ( parmesh->mesh ) {
      fprintf(stdout,"\n  --- MMG ---");
      if ( !MMG2D_defaultValues( parmesh->mesh ) ) {
        ier = 0;
      }
    }
  }

  return ier;
}

int PMMG2D_parsop ( PMMG2D_pParMesh parmesh )
{
  int        ier;

  /* Set mmg verbosity to the max between the Parmmg verbosity and the mmg verbosity */
  assert ( parmesh->mesh->info.imprim == parmesh->info.mmg_imprim );
  parmesh->mesh->info.imprim = MG_MAX ( parmesh->info.imprim, parmesh->mesh->info.imprim );

  parmesh->mesh->info.imprim = parmesh->info.imprim;

  ier = MMG2D_parsop(parmesh->mesh,parmesh->met);

  /* Restore the mesh verbosity */
  parmesh->mesh->info.imprim = parmesh->info.mmg_imprim;

  return ier;
}

void PMMG2D_setfunc( PMMG2D_pParMesh parmesh ) {
  MMG5_pSol met = parmesh->met;

  if( met && met->size == 3 ) {

    PMMG2D_interp3bar = PMMG2D_interp3bar_ani;
    PMMG2D_interp2bar = PMMG2D_interp2bar_ani;

  } else {

    PMMG2D_interp3bar = PMMG2D_interp3bar_iso;
    PMMG2D_interp2bar = PMMG2D_interp2bar_iso;

  }

}
