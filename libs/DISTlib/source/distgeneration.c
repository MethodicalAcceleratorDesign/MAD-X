#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distgeneration.h"
#include "distinterface.h"
#include "outputdist.h"


struct distparam* dist;
struct distparam* diststart;
int dim;
int distn;


void gensixcanonical(){

	int counter = 0;
    int type = dist->incoordtype;
    double tc[dim];
    double normalized[dim], cancoord[dim];
    if(dist->ref->grid==1){
    	generate_grid();//generate a grid
    }

    for(int i =0; i< dist->totincoord; i++){
    	for(int k=0; k<6; k++){
    		tc[k]=dist->incoord[i]->coord[k];

    	}
    	

        if(type==0 || type==3){

        	action2normalized(tc, normalized);
        	normalized2canonical(normalized, cancoord);
        }
        else if(type==1){
	        for(int k=0; k<6; k++){
	    		//normalized[k] = tc[k];
	        	normalized2canonical(tc, cancoord);
	    	}
	    }
	    else if(type==2){
			for(int k=0; k<6; k++){
    			cancoord[k]=tc[k];
    		}

	    }


        if(particle_within_limits_physical(cancoord)==1 && particle_within_limits_normalized(normalized)){
            for(int p=0; p<dim; p++){
                dist->outcoord[counter]->coord[p]   = cancoord[p];

            }
            dist->outcoord[i]->mass  = dist->incoord[i]->mass;
            dist->outcoord[i]->a     = dist->incoord[i]->a;
            dist->outcoord[i]->z     = dist->incoord[i]->z;
            counter++;
       
        }
    }
    dist->totoutcoord=counter;
    dist->isDistrcalculated=1;

}
int gettotalgridlength(){
	int totlength=1;
	for(int i=0; i<dim; i++){
		totlength=totlength*dist->ref->readinlength[i];
	}
	dist->totincoord=totlength;
	return totlength;
}
void generate_grid(){

    int counter = 0;
    double **readin;
    int totallenght = gettotalgridlength();
    readin = (double**)malloc(dim*sizeof(double*));

    for(int i=0; i<dim; i++){
    	readin[i] = (double*)malloc(dist->ref->readinlength[i]*sizeof(double));
    	memcpy(readin[i],dist->incoord[i]->coord, dist->ref->readinlength[i]*sizeof(double));
    }

    deallocateincoord();

    allocateincoord(totallenght);

    //dist->distout_normalized = (double**)malloc(getupperbound()*sizeof(double*));

    for(int i =0; i< dist->ref->readinlength[0]; i++){
        for(int j =0; j< dist->ref->readinlength[1]; j++){
            for(int k =0; k< dist->ref->readinlength[2]; k++){
                for(int l =0; l<dist->ref->readinlength[3]; l++){
                    for(int m =0; m< dist->ref->readinlength[4]; m++){
                        for(int n =0; n< dist->ref->readinlength[5]; n++){
                            dist->incoord[counter]->coord[0]=readin[0][i];
                            dist->incoord[counter]->coord[1]=readin[1][j];
                            dist->incoord[counter]->coord[2]=readin[2][k];
                            dist->incoord[counter]->coord[3]=readin[3][l];
                            dist->incoord[counter]->coord[4]=readin[4][m];
                            dist->incoord[counter]->coord[5]=readin[5][n];
                            //dist->coord[5]->values[n];
                         	counter++;
                        
                                                        
                        }
                    }
                }
            }
        }   
    }
    dist->totincoord=counter;
    dist->isDistrcalculated=1;
}


/*If emittance is defined it converts to canonical coordinates */
void action2normalized(double acangl[6], double normalized[6]){

    normalized[0]= sqrt(acangl[0])*cos(acangl[1]);
    normalized[1]=-sqrt(acangl[0])*sin(acangl[1]);
    normalized[2]= sqrt(acangl[2])*cos(acangl[3]);
    normalized[3]=-sqrt(acangl[2])*sin(acangl[3]);
    normalized[4]= sqrt(acangl[4])*cos(acangl[5]);
    normalized[5]=-sqrt(acangl[4])*sin(acangl[5]); // used to devide with 1000 here before..

}

void normalized2canonical(double normalized_in[6], double cancoord[6]){
    double normalized[6];
    normalized[0] = sqrt(dist->emitt->e1)*normalized_in[0];
    normalized[1] = sqrt(dist->emitt->e1)*normalized_in[1];
    normalized[2] = sqrt(dist->emitt->e2)*normalized_in[2];
    normalized[3] = sqrt(dist->emitt->e2)*normalized_in[3];
    normalized[4] = sqrt(dist->emitt->e3)*normalized_in[4];
    normalized[5] = sqrt(dist->emitt->e3)*normalized_in[5];

       
    if(dist->incoordtype==3) {
        double lindp = 0;
        double lindeltas=0;
        double deltap = normalized[4];
        double deltas = normalized[5];
        double *xap;
        for(int i=0; i<4;i++){
            lindeltas = lindeltas+dist->tas[4][i]*normalized[i];
            lindp=lindp+dist->tas[5][i]*normalized[i];
        } 

        lindp = deltap - lindp;
        lindeltas = deltas - lindeltas;

        xap = (double*)malloc(2*sizeof(double));
        solve2by2eq(dist->tas[4][4], dist->tas[4][5], lindeltas, dist->tas[5][4], dist->tas[5][5], lindp, xap );
        normalized[4] = xap[0];
        normalized[5] = xap[1];

    }
    mtrx_vector_mult_pointer(dim,dim, dist->tas, normalized,cancoord);

}


/*Checks if the particle is within the physical limit set by the user*/
int particle_within_limits_physical(double *physical){
    
    if(dist->cuts2apply->isset_p==0) return 1;
    for(int i=0; i<dim; i++){
        if(dist->cuts2apply->physical[i]->isset==1){
            if(physical[i] > dist->cuts2apply->physical[i]->min && physical[i] < dist->cuts2apply->physical[i]->max) return 0;
        }   
    }
    
    return 1;

}

int particle_with_limits_action(int i, double value){
    
    if(dist->cuts2apply->action[i]->isset==1){
        if(value > pow(dist->cuts2apply->action[i]->min,2) && value < pow(dist->cuts2apply->action[i]->max,2)) return 1;
        else return 0;
    }
    else{
        return 1;
    }

    return 1;
}

/*Checks if the particle is within the normalized limit set by the user*/
int particle_within_limits_normalized(double *normalized){
    
    if(dist->cuts2apply->isset_n==0) return 1;
    for(int i=0; i<dim; i++){
        if(dist->cuts2apply->normalized[i]->isset==1){
            if(normalized[i] > dist->cuts2apply->normalized[i]->min && normalized[i] < dist->cuts2apply->normalized[i]->max) return 0;
        }
    }
    return 1;
}

void createcoordinates(int index,  double start, double stop, int length, int type){
	double *temp;
    temp = (double*)malloc(length*sizeof(double));
    
    if(type ==0){ //Constant value
    	for(int i=0;i <length; i++){
        	dist->incoord[i]->coord[index] = start;
    	}
    return;
    }
   

    else if(type==1){ //Linearly spaced intervalls
		createLinearSpaced(length, start, stop, temp);
        for(int i=0;i <length; i++){
            dist->incoord[i]->coord[index] = temp[i];
        }

    }
    else if(type==2){ //Exponentially spaced
    	createLinearSpaced(length, start, stop, temp);
        for(int i=0;i <length; i++){
            dist->incoord[i]->coord[index] = exp(temp[i]);
        }
    }
    else if(type==3){ //Spaced with  ^2
    	createLinearSpaced(length, start, stop, temp);
        for(int i=0;i <length; i++){
            dist->incoord[i]->coord[index] = pow(temp[i],2);
        
        }

    }
    else if(type==4){ // uniform random 
        for(int i=0;i <length; i++){
            dist->incoord[i]->coord[index] = rand_uni(start, stop);
        }
    }

    else if(type==5){ // Gaussian random (Here start is mean and stop is the standard deviation)
        for(int i=0;i <length; i++){
            dist->incoord[i]->coord[index] = randn(start, stop);
        }
    }

    else if(type==6){ // Rayleigh distribution
        double tmp;
        for(int i=0;i <length; i++){
            
            tmp = randray(start, stop);
            while(particle_with_limits_action(index, tmp)==0)
                tmp = randray(start, stop);

            dist->incoord[i]->coord[index] = tmp;
        }
    }
    else
    	issue_error("Unknown type of spacing");

    free(temp);
}


void createLinearSpaced(int length, double start, double stop, double *eqspaced ){
    
    double distance = (stop-start)/length;
    for(int i=0; i<length; i++){
        eqspaced[i] = start+distance*i;
    }
}

double randn(double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ( rand () / (double)RAND_MAX) * 2;
      U2 = -1 + ( rand () / (double)RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

double rand_uni(double low, double high)
{

  return (double)(rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

double randray(double mu, double sigma){
  double low = 0;
  double high =1;
  return pow((mu+sigma*sqrt((-2*log(rand_uni(low, high))))),2);
}


void allocateincoord(int linecount){
  dist->incoord = (struct coordinates**)malloc(linecount*sizeof(struct coordinates*));
  dist->outcoord = (struct coordinates**)malloc(linecount*sizeof(struct coordinates*));
  dist->totincoord = linecount;
  for(int i=0; i< dim; i++){
    dist->ref->typeused[i]=-1;
  }
  for(int i=0; i<linecount; i++){
    dist->incoord[i] = (struct coordinates*)malloc(sizeof(struct coordinates));
    
    dist->incoord[i]->coord = (double*)malloc(dim*sizeof(double));
    dist->incoord[i]->readin = (double*)malloc(32*sizeof(double));

    dist->incoord[i]->mass=0;
    dist->incoord[i]->a=0;
    dist->incoord[i]->z=0;

    dist->outcoord[i] = (struct coordinates*)malloc(sizeof(struct coordinates));
    dist->outcoord[i]->coord = (double*)malloc(dim*sizeof(double));

  }
  dist->isallocated =1;
}
void deallocateincoord(void){


  free(dist->incoord); 
  free(dist->outcoord); 
  //dist->ref->typeused = (int*)malloc(dim*sizeof(int));
  dist->isallocated = 0;
  dist->totincoord = -1;
}