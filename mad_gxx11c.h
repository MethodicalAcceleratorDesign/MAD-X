#ifndef MAD_GXX11C_H
#define MAD_GXX11C_H

// Called by gxx11.f90
void wopen   (int *uswid, int *ushi);
void wclose  (void);
void wclrwk  (int *i1,int *i2);
void wpl     (int *np, float *xp, float *yp);
void wfa     (int *np, float *xp, float *yp);
void wswn    (float *wlx, float *wxfact, float *wly, float *wyfact);
void wtx     (float *xp,float *yp, char *string);
void wwait   (void);
void wsetci  (char *uscol);
void wsetls  (int *ls);
void wstring (char *s, int *l);
void cbyt    (int* source, int* s_pos, int* target, int* t_pos, int* n);
void mydtime (int* year, int* month, int* day, int* hour, int* minute, int* sec);

#endif // MAD_GXX11C_H

