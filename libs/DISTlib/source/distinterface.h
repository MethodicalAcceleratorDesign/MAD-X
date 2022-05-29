struct cut{
	int isset;
	double min;
	double max;
};


void writefile_f(const char*  filename_in, int strlen);
void readfile_f(const char*  filename_in, int strlen);
void readtasmatrixfile(const char*  filename_in);
void initializedistribution(int numberOfDist);
void setdistribution(int ndist);
void setphysicalcut(int variable, double min, double max);
void setnormalizedcut(int variable, double min, double max);

void settotalsteps(int totgenerate);
void sete0andmass0(double energy0, double mass0);
void setemitt12(double e1, double e2);
void setemitt3(double e3);
void settasmatrix(double * tas);
void settasmatrixtranspose(double *tas);
void setcoords(double *xn, double *xnp, double *yn, double *ynp, double *zn, double *znp, int totparticles, int coordtype);
void get6trackcoord(double *x, double *xp, double *y, double *yp, double *sigma, double *deltap, int *totparticles);
void setscan_para_diagonal(int variable, int variable_type, int type, double start, double stop);
void setscan_para_grid(int variable,int variable_type,int type, double start, double stop, int length);
void setdisptas(double dx, double dpx, double dy, double dpy);
void settwisstas(double betax, double alfax, double betay, double alfay);
void getunconvertedcoord(double *x, double *xp, double *y, double *yp, double *sigma, double *deltap, int *totparticles);
void setactionanglecut(int variable, double min, double max);
void free_distribution(int i);
void settasmatrix_element(double value, int row, int column);
void addclosedorbit(double *clo);
void getarraylength(int *totlength);
void getrefpara(double *energy0, double *mass0, int *a0, int *z0);
