#define MAX_COLUMNS 100
#define MAX_LENGTH 100
#define MAX_ROWS 100000
void canonical2six(double *canonical, double beta0, double pc0, double mass, double *coord);
double momentum2energy(double momentum, double mass);
double energy2momentum(double energy, double mass);
double pt2deltap(double pt, double beta0 );
double psigma2deltap(double psigma, double beta0 );
void issue_error(const char* t1);
void issue_error2(const char* t1,const char* t2);
void issue_warning(const char* t1);
void issue_info(const char* t1);
int strcmpnl (const char *s1, const char *s2);
double sigma2zeta(double sigma, double beta0, double beta);
double tau2zeta(double tau, double beta);
void solve2by2eq(double a1, double b1, double c1, double a2, double b2, double c2, double *x);
void mtrx_vector_mult_pointer(int mp, int np,  double **mtrx_a, double mtrx_b[6], double result[6]);
void checkiftimeset(int entype);
void checkifangset(int entype);
void checkifenergyset(int entype);