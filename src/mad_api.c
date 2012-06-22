#include "madx.h"
#include "mad_api.h"

struct sequence_list*
madextern_get_sequence_list(void)
{
    return sequences;
}
