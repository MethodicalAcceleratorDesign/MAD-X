#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinterface.h"
#include "distgeneration.h"
#include "outputdist.h"
#include "distdata.h"


/*
This allocates the the memory for the distributions
*/
void initializedistribution(int numberOfDist){
    dist = (struct distparam*)malloc((numberOfDist)*sizeof(struct distparam));
    dim  = 6;
    
        for(int i = 0; i <numberOfDist; i++)
        {
        
        (dist + i)->ref = (struct refparam*)malloc(sizeof(struct refparam));
        (dist + i)->coord = (struct parameters**)malloc(dim*sizeof(struct parameters*));
        (dist + i)->emitt = (struct emittances*)malloc(sizeof(struct emittances));
        (dist + i)->cuts2apply = (struct appliedcut*)malloc(sizeof(struct appliedcut));
        (dist + i)->cuts2apply->physical = (struct cut**)malloc(dim*sizeof(struct cut*));
        (dist + i)->cuts2apply->normalized = (struct cut**)malloc(dim*sizeof(struct cut*));
        (dist + i)->cuts2apply->action = (struct cut**)malloc(dim*sizeof(struct cut*));
        (dist + i)->tas   = (double**)malloc(dim*sizeof(double*));
        (dist + i)->invtas   = (double**)malloc(dim*sizeof(double*));
        (dist + i)->closedorbit   = (double*)malloc(dim*sizeof(double));
        (dist + i)->isDistrcalculated = 0;
		(dist + i)->ref->e0=0;
		(dist + i)->ref->pc0=0;
		(dist + i)->ref->a0=1;
		(dist + i)->ref->z0=1;
		(dist + i)->ref->mass0=0;
		(dist + i)->ref->charge0=1;
		(dist + i)->ref->en_like=-1;
        (dist + i)->ref->time_like=-1;
        (dist + i)->ref->ang_like=-1;
        (dist + i)->totincoord=-1;
        (dist + i)->ref->typeused = (int*)malloc(dim*sizeof(int));
        (dist + i)->ref->readinlength = (int*)malloc(dim*sizeof(int));
        (dist + i)->ref->grid=-1;
        for(int k=0; k<dim;k++){
            (dist + i)->tas[k] =(double*)malloc(dim*sizeof(double));
            (dist + i)->invtas[k] =(double*)malloc(dim*sizeof(double));
        }
        (dist + i)->incoordtype   =-1;
        (dist + i)->disttype = 0;
        for(int j=0; j<dim; j++)
        {
            (dist + i)->cuts2apply->physical[j]   = (struct cut*)malloc(sizeof(struct cut));
            (dist + i)->cuts2apply->normalized[j] = (struct cut*)malloc(sizeof(struct cut));
            (dist + i)->cuts2apply->action[j]     = (struct cut*)malloc(sizeof(struct cut));
            (dist +i)->coord[j] = (struct parameters*)malloc(sizeof(struct parameters));
            (dist +i)->coord[j]->start=0;
            (dist +i)->coord[j]->stop=0;
            (dist +i)->coord[j]->length=1;
            (dist +i)->coord[j]->type=0;
            (dist +i)->closedorbit[j]=0;
        }
        for (int i = 0; i < 6; i++)
        {
            for (int k = 0; k < 6; k++)
            {
                dist->tas [i][k] = 0;
                dist->invtas[i][k] = 0;
                
            }
        }
    }
    diststart=dist;

}

void sete0andmass0(double energy0, double mass0){
	dist->ref->mass0 = mass0;
	dist->ref->e0 = energy0;
    dist->ref->pc0 = energy2momentum(dist->ref->e0,dist->ref->mass0);
    dist->ref->beta0 = (dist->ref->pc0)/(dist->ref->e0);
}

void setdistribution(int ndist){
		dist = diststart + ndist;
}

void setemitt12(double e1, double e2){
    dist->emitt->e1=e1; 
    dist->emitt->e2=e2; 
}

void setemitt3(double e3){
        dist->emitt->e3=e3; 
}

void settasmatrix(double *tas){
	for(int i =0; i<dim; i++){
		for(int j =0; j<dim; j++){
			dist->tas[i][j] = tas[j+i*dim];
		}
	}
}

void settasmatrixtranspose(double *tas){
    for(int i =0; i<dim; i++){
        for(int j =0; j<dim; j++){
            dist->tas[j][i] = tas[j+i*dim];
        }
    }
}

void setactionanglecut(int variable, double min, double max){
    dist->cuts2apply->isset_a=1;
    dist->cuts2apply->action[variable]->min=min;
    dist->cuts2apply->action[variable]->max=max;
    dist->cuts2apply->action[variable]->isset=1;
}

void settasmatrix_element(double value, int row, int column){

            dist->tas[row][column] = value;
        
}
/*Create a TAS matrix using E-T, assuming uncoupled system */
void settwisstas(double betax, double alfax, double betay, double alfay){

    dist->tas[0][0] = sqrt(betax);
    dist->tas[1][0] =-alfax/sqrt(betax);
    dist->tas[1][1] =-1/sqrt(betax);

    
    dist->tas[2][2] = sqrt(betay);
    dist->tas[3][2] =-alfay/sqrt(betay);
    dist->tas[3][3] =-1/sqrt(betay);

    dist->tas[5][5] =1; // This is not a matched longitudinal but enables to change the energy 
    dist->tas[4][4] =1;

}

void setdisptas(double dx, double dpx, double dy, double dpy){
    
    dist->tas[0][5] = dx;
    dist->tas[1][5] = dpx;
    dist->tas[2][5] = dy;
    dist->tas[3][5] = dpy;

}


// 0 -action angle, 1- normalized coordinates, 2-physical, 3-mixed
void setcoords(double *xn, double *xpn, double *yn, double *ypn, double *zn, double *zpn, int totparticles, int coordtype){
	allocateincoord(totparticles);
	for(int i=0; i<totparticles; i++){
		dist->incoord[i]->coord[0] = xn[i];
		dist->incoord[i]->coord[1] = xpn[i];
		dist->incoord[i]->coord[2] = yn[i];
		dist->incoord[i]->coord[3] = ypn[i];
		dist->incoord[i]->coord[4] = zn[i];
		dist->incoord[i]->coord[5] = zpn[i]; 
    }

	dist->totincoord  = totparticles;
	dist->incoordtype = coordtype;
}

void settotalsteps(int totgenerate){
    dist->totincoord = totgenerate;
}

void setscan_para_diagonal(int variable, int variable_type, int type, double start, double stop){
    if(dist->totincoord==-1)
        issue_error("Must set a total steps before you can set a parameter");

    if(dist->isallocated!=1){
        allocateincoord(dist->totincoord);
    }
    dist->ref->typeused[variable] =  variable_type;
    dist->incoordtype=variable_type; // this is wrong
    
    createcoordinates(variable, start, stop, dist->totincoord,type);
}

void setscan_para_grid(int variable,int variable_type,int type, double start, double stop, int length){
    dist->ref->readinlength[variable] = length;

    if(length>MAX_LENGTH)
        issue_error("For grid scans you have requested to many particles.");
    if(dist->isallocated!=1)
        allocateincoord(MAX_LENGTH);

    dist->ref->typeused[variable] =  variable_type;
    dist->incoordtype=variable_type;
    createcoordinates(variable, start, stop, dist->totincoord,type);
    dist->ref->grid=1;

}

void addclosedorbit(double *clo){
	for(int i=0; i<dim;i++){
		dist->closedorbit[i] = clo[i];
    }
}


void setphysicalcut(int variable, double min, double max){
	dist->cuts2apply->isset_p=1;
	dist->cuts2apply->physical[variable-1]->min=min;
	dist->cuts2apply->physical[variable-1]->max=max;
	dist->cuts2apply->physical[variable-1]->isset=1;

}

void setnormalizedcut(int variable, double min, double max){
	dist->cuts2apply->isset_n=1;
	dist->cuts2apply->normalized[variable-1]->min=min;
	dist->cuts2apply->normalized[variable-1]->max=max;
	dist->cuts2apply->normalized[variable-1]->isset=1;

}
void getarraylength(int *totlength){
    if(dist->isDistrcalculated ==0){

        gensixcanonical();
    }
    *totlength=dist->totoutcoord;
}

void get6trackcoord(double *x, double *xp, double *y, double *yp, double *sigma, double *deltap, int *totparticles){
    double tmp[6];
    int nparticles;
    if(dist->isDistrcalculated ==0){
        gensixcanonical();

    }
    if(dist->totoutcoord < *totparticles)
        nparticles = dist->totoutcoord;
    else
        nparticles = *totparticles;

    for(int i=0; i < nparticles; i++){
        canonical2six(dist->outcoord[i]->coord, dist->ref->beta0, dist->ref->pc0, dist->incoord[i]->mass, tmp);
        x[i]  = tmp[0];
        xp[i] = tmp[1];
        y[i]  = tmp[2];
        yp[i] = tmp[3];
        sigma[i]  = tmp[4];
        deltap[i] = tmp[5];
   }
    *totparticles=nparticles;
}

void getunconvertedcoord(double *x, double *xp, double *y, double *yp, double *sigma, double *deltap, int *totparticles){
    int nparticles;
    if(dist->isDistrcalculated ==0){
        gensixcanonical();

    }
    if(dist->totoutcoord < *totparticles)
        nparticles = dist->totoutcoord;
    else
        nparticles = *totparticles;

    for(int i=0; i < nparticles; i++){
        x[i]  = dist->outcoord[i]->coord[0];
        xp[i] = dist->outcoord[i]->coord[1];
        y[i]  = dist->outcoord[i]->coord[2];
        yp[i] = dist->outcoord[i]->coord[3];
        sigma[i]  = dist->outcoord[i]->coord[4];
        deltap[i] = dist->outcoord[i]->coord[5];
   }
    *totparticles=nparticles;
}

void readtasmatrixfile(const char*  filename_in){
    FILE * fp;
    char line [500]; /* or other suitable maximum line size */

    fp = fopen(filename_in, "r+");
    int i=0;
    if(fp==NULL){
        printf("No such file\n");
        exit(1);
    }
   while ( fgets ( line, sizeof line, fp ) != NULL ){ /* read a line */
        printf("%s \n", line);
        sscanf(line, "%lf %lf %lf %lf %lf %lf", &dist->tas[i][1], &dist->tas[i][1], &dist->tas[i][2], &dist->tas[i][3], &dist->tas[i][4], &dist->tas[i][5]);
        //printf("%f %f %f %f %f %f \n", dist->tas[i][1], dist->tas[i][1], dist->tas[i][2], dist->tas[i][3], dist->tas[i][4], dist->tas[i][5]);
        i++;
    }
    fclose(fp);
}


void getrefpara(double *energy0, double *mass0, int *a0, int *z0){
    *energy0= dist->ref->e0;
    *mass0  = dist->ref->mass0;
    *a0     = dist->ref->a0;
    *z0     = dist->ref->z0;
}
/*
void readfile_f(const char*  filename_in, int strlen){
    char filename [strlen];
    strncpy(filename, filename_in, strlen);
    filename[strlen] = '\0';
    readfile(filename);
}*/

void writefile_f(const char*  filename_in, int strlen){
    char filename [strlen];
    strncpy(filename, filename_in, strlen);
    filename[strlen] = '\0';
    print2file(filename);
}

void free_distribution(int i){

    free((dist+i)->closedorbit);
    free((dist+i)->tas);
    free((dist+i)->invtas);

    for(int j=0; j<dim; j++)
    {
        free((dist + i)->cuts2apply->physical[j]) ;
        free((dist + i)->cuts2apply->normalized[j]) ;
        free((dist + i)->cuts2apply->action[j]) ;
        free((dist +i)->coord[j]);
    }
    free((dist+i)->coord);
    free((dist + i)->cuts2apply->physical);
    free((dist + i)->cuts2apply->normalized);
    free((dist + i)->cuts2apply->action);
    free((dist + i)->cuts2apply);
    
    free((dist + i)->ref);
    free((dist + i)->emitt);
    free(dist+i);


}