#include "madx.h"
#include "mad_api.h"

void
madextern_start(void)
{
    madx_start();
}

struct sequence_list*
madextern_get_sequence_list(void)
{
    return sequences;
}
void madextern_end()
{
    madx_finish();
}
void madextern_input(char* ch) 
{
    stolower_nq(ch);
    pro_input(ch);
}

