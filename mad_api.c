#include "mad_api.h"
#include "mad_str.h"
#include "mad_core.h"
#include "mad_eval.h"
#include "mad_wrap_f.h"

// dependencies from madxp.c

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

