#include "mad_math.h"

double
fact(int i)
{
  int k;
  double nfact;

  if(i == 0) {
     return(1.0);
  } else if(i == 1) {
     return(1.0);
  } else if(i > 1) {
    nfact = 1.;
    for (k=1;k<=i;k++) {
        nfact = nfact*k;
    }
     return(nfact);
  } else if(i < 0){
     return(-1.0);
  } else {
     return(-1.0);
  }
}


