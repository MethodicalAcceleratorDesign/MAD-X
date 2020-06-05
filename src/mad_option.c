#include "madx.h"

 int
get_option(const char* str)
{
/* This function is called by fortran to get option of a command */
  int i;
  mycpy(c_dum->c, str);
  if (options != NULL && (i = name_list_pos(c_dum->c, options->par_names)) > -1)
    return options->par->parameters[i]->double_value; // (k = not used
  else if (strcmp(c_dum->c, "warn") == 0) return init_warn;
  else return 0;
}

void
set_option(const char* str, int* opt)
  /* sets an (old or new) option with name "str",
     value *opt (0 false, 1 true) */
{
  int i, k;
  char* bc;
  mycpy(c_dum->c, str); bc = permbuff(c_dum->c);
  if ((i = name_list_pos(bc, options->par_names)) < 0)
  {
    add_to_name_list(bc, 0, options->par_names); // j = not used
    if ((k = options->par->curr) == options->par->max)
      grow_command_parameter_list(options->par);
    options->par->parameters[options->par->curr++] = new_command_parameter(bc, 0);
    options->par->parameters[k]->double_value = *opt;
  }
  else options->par->parameters[i]->double_value = *opt;
}

void
set_defaults(const char* str) /* reset options, beam etc. to defaults */
{
  int i, pos;
  struct command* beam_clone;

  if ((pos = name_list_pos(str, defined_commands->list)) > -1)
  {
    if (strcmp(str, "option") == 0)
    {  
      if (options != NULL) delete_command(options);
      options = clone_command(defined_commands->commands[pos]);
    }
    else if (strcmp(str, "set") == 0)
      store_set(defined_commands->commands[pos], 0);
    else if (strcmp(str, "setplot") == 0)
    {
      if (plot_options != NULL) delete_command(plot_options);
      plot_options = clone_command(defined_commands->commands[pos]);
    }
    else if (strcmp(str, "threader") == 0)
    {
      if (threader_par != NULL)  delete_command(threader_par);
      threader_par = clone_command(defined_commands->commands[pos]);
    }
    else if (strcmp(str, "beam") == 0)
    {
      if (current_beam == NULL)
        current_beam = clone_command(defined_commands->commands[pos]);
      beam_clone = clone_command(defined_commands->commands[pos]);
      for (i = 0; i < beam_clone->par_names->curr; i++)
        beam_clone->par_names->inform[i] = 1; /* mark as "read" */
      update_beam(beam_clone);
      delete_command(beam_clone);
    }
  }
}

