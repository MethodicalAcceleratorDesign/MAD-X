#ifndef MAD_RAND_H
#define MAD_RAND_H

void    init55(int seed);
double  frndm(void);
double  grndm(void);
double  tgrndm(double cut);

void    setrand (const char *kind, int rng_id); // "default" or "best"

#endif // MAD_RAND_H


