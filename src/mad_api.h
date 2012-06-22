#ifndef MAD_API_H
#define MAD_API_H

/**
 * @file mad_api.h Calls to use for external libraries
 *
 * This is a header that can be included in external libraries that wants to 
 * link against the mad-x library. It was mostly needed because
 * we did not have clean headers before. 
 * 
 * Since that has now changed, have a look instead at the functions 
 * madx_start() and madx_finish() in mad_core.h, 
 * and pro_input() in mad_eval.h. Those provide everything
 * you need for a basic interface.
 * 
 * @author Yngve Inntjore Levinsen
 */

// types

struct sequence_list;

// variables

extern struct sequence_list* sequences;

/**
 * @brief Get list of available sequences
 *
 * This function is part of the madextern external functions.
 * These are not supposed to be used by any internal mad-x code.
 *
 * This function returns a name_list structure of all currently
 *  available sequences. Reason we do not return sequence_list
 *  is that name_list are "cleaner" structures which are easier
 *  to access in python (currently only use-case).
 *
 */
struct sequence_list* madextern_get_sequence_list(void);

#endif // MAD_API_H

