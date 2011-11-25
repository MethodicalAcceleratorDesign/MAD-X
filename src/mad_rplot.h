#ifndef MAD_RPLOT_H
#define MAD_RPLOT_H

/* Piotr Skowronski, CERN */

void type_ofCall newrplot(void);
void type_ofCall plottrack(int* particleno, int* obspoint, int* turn,
                           double* x, double* xp,
                           double* y, double* yp,
                           double* dpOverP,  double* p,
                           double* length);
void type_ofCall plottwiss(int* obspoint,
                           double* betax, double* alfax,
                           double* betay, double* alfay,
                           double* betaz, double* alfaz,
                           double* length);
void type_ofCall rplotfinish(void);
void type_ofCall rviewer(void);
void type_ofCall madxv_setfctnname  (int* n, const char* name);
void type_ofCall madxv_setknobname  (int* n, const char* name);
void type_ofCall madxv_setfunctionat(int* n, int* el, const char* name);

#endif // MAD_RPLOT_H
