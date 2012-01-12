#ifndef MAD_RPLOT_H
#define MAD_RPLOT_H

/* Piotr Skowronski, CERN */

void newrplot(void);
void plottrack(int* particleno, int* obspoint, int* turn,
                    double* x, double* xp,
                    double* y, double* yp,
                    double* dpOverP,  double* p,
                    double* length);
void plottwiss(int* obspoint,
               double* betax, double* alfax,
               double* betay, double* alfay,
               double* betaz, double* alfaz,
               double* length);
void rplotfinish(void);
void rviewer(void);
void madxv_setfctnname  (int* n, const char* name);
void madxv_setknobname  (int* n, const char* name);
void madxv_setfunctionat(int* n, int* el, const char* name);

#endif // MAD_RPLOT_H
