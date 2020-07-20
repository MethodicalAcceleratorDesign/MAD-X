#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinterface.h"
#include "distgeneration.h"
#include "outputdist.h"

#include "distdata.h"


void print2file(const char* nameoffile){

   FILE * fp;
   /* open the file for writing*/
   fp = fopen (nameoffile,"w");
   fprintf (fp, "@ mass0 %f \n",dist->ref->mass0);
   fprintf (fp, "@ charge0 %d \n",dist->ref->charge0);
   fprintf (fp, "@ z0 %d \n",dist->ref->z0);
   fprintf (fp, "@ a0 %d \n",dist->ref->a0);
   fprintf (fp, "@ pc0 %f \n",dist->ref->pc0);

   fprintf(fp, "x px y py zeta deltap \n");
   if(dist->isDistrcalculated ==0){
        gensixcanonical();
    }
   for(int i=0; i<dist->totoutcoord; i++){
   	fprintf(fp, "%.9e %.9e %.9e %.9e %.9e %.9e \n", dist->outcoord[i]->coord[0],dist->outcoord[i]->coord[1],dist->outcoord[i]->coord[2],
   		dist->outcoord[i]->coord[3],dist->outcoord[i]->coord[4],dist->outcoord[i]->coord[5]);
   }

   /* close the file*/  
   fclose (fp);
}