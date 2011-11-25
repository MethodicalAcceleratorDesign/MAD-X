#include "madx.h"

void
pro_survey(struct in_cmd* cmd)
  /* calls survey module */
{
  struct name_list* nl = current_survey->par_names;
  struct command_parameter_list* pl = current_survey->par;
  char *filename = NULL, *table_name;
  int pos, w_file;
  int iarc = 1, keep;

  (void)cmd;
  if (current_sequ == NULL)
  {
    warning("SURVEY, but no active sequence:", "ignored");
    return;
  }
  if (debuglevel > 1) fprintf(prt_file, "enter Survey module\n");
  keep = get_option("rbarc");
  set_option("rbarc", &iarc);
  pos = name_list_pos("file", nl);
  if (nl->inform[pos])
  {
    if ((filename = pl->parameters[pos]->string) == NULL)
    {
      if (pl->parameters[pos]->call_def != NULL)
        filename = pl->parameters[pos]->call_def->string;
    }
    if (filename == NULL) filename = permbuff("dummy");
    w_file = 1;
  }
  else w_file = 0;
  pos = name_list_pos("table", nl);
  if(nl->inform[pos]) /* table name specified - overrides save */
  {
    if ((table_name = pl->parameters[pos]->string) == NULL)
      table_name = pl->parameters[pos]->call_def->string;
  }
  else table_name = permbuff("survey");
  survey_table = make_table(table_name, "survey", survey_table_cols,
                            survey_table_types, current_sequ->n_nodes);
  add_to_table_list(survey_table, table_register);
  survey_();
  if (w_file) out_table(table_name, survey_table, filename);
  set_option("rbarc", &keep);
}


