#include "madx.h"

void
pro_match(struct in_cmd* cmd)
  /* controls the matching module */
{
  /* OB 12.2.2002: changed the sequence of if statements so that MAD
     can go through the whole matching sequence */

  if (strcmp(cmd->tok_list->p[0], "match") == 0)
  {
    match_match(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "cell") == 0)
  {
    warning("CELL command no longer valid, ","use MATCH");
    return;
  }
  else if (match_is_on == 0)
  {
    warning("no MATCH command seen,","ignored");
    return;
  }
  else if (strcmp(cmd->tok_list->p[0], "endmatch") == 0)
  {
    match_end(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "migrad") == 0 ||
           strcmp(cmd->tok_list->p[0], "lmdif") == 0 ||
           strcmp(cmd->tok_list->p[0], "jacobian") == 0 ||
           strcmp(cmd->tok_list->p[0], "simplex") == 0)
  {
    match_action(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "constraint") == 0)
  {
    match_constraint(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "couple") == 0)
  {
    match_couple(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "fix") == 0)
  {
    match_fix(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "global") == 0)
  {
    match_global(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "level") == 0)
  {
    match_level(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "vary") == 0)
  {
    match_vary(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "weight") == 0)
  {
    match_weight(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "gweight") == 0)
  {
    match_gweight(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "rmatrix") == 0)
  {
    match_rmatrix(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "tmatrix") == 0)
  {
    match_tmatrix(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "global") == 0)
  {
    match_global(cmd);
  }
  else if (strcmp(cmd->tok_list->p[0], "use_macro") == 0)
  {
    match2_macro(cmd);
  }
}

