struct distparam
{
	struct parameters** coord;
    int disttype;
	struct refparam* ref;
	struct emittances* emitt;
	double **tas;
	double **invtas;
	double *closedorbit;
	int incoordtype; // This tells which type of coordinates the input is given.  // 0-action angle, 1-Normalized, 2-physical, 3-mixed action angle and physical (lat)
	struct coordinates** incoord;
	struct coordinates** outcoord;
	struct coordinates** gridin;
	int totincoord;
	int totoutcoord;
	int isDistrcalculated;
	int isallocated;
	struct appliedcut* cuts2apply;
};
struct refparam{
	int charge0;
	int z0;
	int a0;
	double pc0;
	double e0;
	double mass0;
	double beta0;
	int en_like;
	int time_like;
	int ang_like;
	struct emittances* emitt;
	int *typeused; // This says which is used or a mixture of it.. (used for reading in and as a cross check)
	// 0-action, 1-normalized, 2-physical
	int *readinlength;
	int grid;

};
struct coordinates
{
	double *coord;
	double *readin;
	double mass;
	int charge;
	int a;
	int z;

};

struct parameters
{
  double start;
  double stop;   
  int length;
  int type; //This gives the type of distribution, constant, linear, gaussian, 
  double * values;
  int coordtype; 
};

struct emittances{
	double e1, e2, e3;
	double dp, deltas;
};

struct appliedcut{
	int isset_p;
	int isset_n;
	int isset_a;
	struct cut** physical;
	struct cut** normalized;
	struct cut** action;
};

//void canonical2emittance_(double cancord[6], double emittance[3]);
//void dist2sixcoord_();
//void action2canonical_(double acangl[6], double cancord[6], double acoord[6]);
//void action2sixinternal_(double tc[6], double *results, double *normalized);
//int checkdist();
//void createrandomdist_();
void gensixcanonical(void);
int particle_within_limits_physical(double *physical);
int particle_within_limits_normalized(double *normalized);
int particle_with_limits_action(int i, double value);
void generatefromnormalized(void);
void generatefrommixed(void);
void generatefromaction(void);
void generatefromphysical(void);
void action2normalized(double acangl[6], double normalized[6]);
void normalized2canonical(double normalized[6], double cancoord[6]);
double randn(double mu, double sigma);
double rand_uni(double low, double high);
void createLinearSpaced(int length, double start, double stop, double *eqspaced );
double randray(double mu, double sigma);
void createcoordinates(int index,  double start, double stop, int length, int type);
int gettotalgridlength(void);
void generate_grid(void);
void allocateincoord(int linecount);
void deallocateincoord(void);