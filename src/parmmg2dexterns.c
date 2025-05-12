#include "parmmg2dexterns.h"
#include "parmmg2d.h"

int (*PMMG2D_interp4bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTetra pt,int,PMMG2D_barycoord*)=NULL;
int (*PMMG2D_interp3bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int,PMMG2D_barycoord*)=NULL;
int (*PMMG2D_interp2bar)(MMG5_pMesh mesh,MMG5_pSol met,MMG5_pSol oldMet,MMG5_pTria ptr,int ip,int l,PMMG2D_barycoord *barycoord)=NULL;
