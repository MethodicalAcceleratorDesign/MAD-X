#include <unistd.h>
#include <stdio.h>
#include <signal.h>

void doagain(int dummy)
{
}

void epause_()
{
        int i;

/*      fprintf(stderr," PAUSE called! \n"); */
        signal(SIGCONT,*doagain);
        i=pause();
/*      fprintf(stderr," PAUSE returned! \n"); */
}
