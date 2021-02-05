#include "madx.h"

static void
get_pos_slice_in_node(struct node* cnode, double pos, double* pos_vector){

  struct node* tmpnode;
  tmpnode = current_node;
  current_node = cnode;
  printf("iii %f \n", node_value("dx"));
  locslice_(&pos, pos_vector);

  current_node = tmpnode;

}


static void
test_align_command(void){
  //This was just used for debug and will be fixed later... 
  struct node* node = current_sequ->start;
  double tmp[7];
  while (node != NULL) {
    get_pos_slice_in_node(node, 0.0, tmp);
    printf("%s x: %f y: %f s: %f \n", node->name, tmp[0], tmp[1], tmp[2] );
    get_pos_slice_in_node(node, 0.1, tmp);
    printf("%s x: %f y: %f s: %f \n", node->name, tmp[0], tmp[1], tmp[2] );
    get_pos_slice_in_node(node, 0.2, tmp);
    current_node = node;
    printf("%s x: %f y: %f s: %f \n", node->name, tmp[0], tmp[1], tmp[2] );
    if (node == current_sequ->end) break;
    node = node->next;
  }
}

void
pro_survey(struct in_cmd* cmd)
  /* calls survey module */
{
  struct sequence* keep_current;
  char *filename = NULL, *table_name;
  int w_file;
//  int iarc = 1, keep;

  (void)cmd;
  if (current_sequ == NULL)
  {
    warning("SURVEY, but no active sequence:", "ignored");
    return;
  }
  if (debuglevel > 1) fprintf(prt_file, "enter Survey module\n");
//  LD: why is Survey blocking the RBARC option?
//  RBARC is used by el_par_value in mad_elem.c

//  keep = get_option("rbarc");
//  set_option("rbarc", &iarc);
  w_file = command_par_string_user2("file", current_survey, &filename);
  if (w_file && !filename)
    filename = permbuff("dummy");
  table_name = command_par_string_user("table", current_survey);
  if(!table_name) table_name = permbuff("survey");
  survey_table = make_table(table_name, "survey", survey_table_cols,
                            survey_table_types, current_sequ->n_nodes);
  add_to_table_list(survey_table, table_register);
  keep_current = current_sequ;
  int l_survey = command_par_value("perm_align_survey", current_survey);
  
  if(l_survey==1) elementloc_();
  else            survey_();

  current_sequ = keep_current;
  if (w_file) out_table(table_name, survey_table, filename);
// set_option("rbarc", &keep);
  if(debuglevel >1 ) test_align_command();
}

void
pro_use_survey(void)
{
  /* Constructs artificial survey command for USE,SURVEY.
     The survey data are stored at the nodes. */
  /* 2013-Jul-18  19:17:06  ghislain: DOC undocumented feature ? */
  struct in_cmd* pro_use = new_in_cmd(10);
  struct name_list* usenl;
  int usepos;
  pro_use->label = NULL;
  pro_use->type = 0;
  pro_use->clone = pro_use->cmd_def = clone_command(find_command("survey",defined_commands));
  usenl = pro_use->cmd_def->par_names;
  usepos = name_list_pos("table", usenl);
  pro_use->cmd_def->par->parameters[usepos]->string = tmpbuff("survey");
  pro_use->cmd_def->par_names->inform[usepos] = 1;
  usepos = name_list_pos("file", usenl);
  pro_use->cmd_def->par->parameters[usepos]->string = NULL;
  pro_use->cmd_def->par_names->inform[usepos] = 0;
  current_survey=(pro_use->clone);
  pro_survey(pro_use);
  exec_delete_table("survey");





}

