#include "madx.h"

static void
move_files(const char* orig_name,const  char* append,const  char* dirname){
  char dest_name[100];
  strcpy(dest_name, dirname);
  strcat(dest_name, "/");
  strcat(dest_name, orig_name);
  strcat(dest_name, append);
  rename(orig_name, dest_name);
}

static void copy_input_script(const char* filename){
  FILE *source, *target;
  char ch;

  source = fopen(mad_argv[1], "r");
  target = fopen(filename, "w");

  while ((ch = fgetc(source)) != EOF)
      fputc(ch, target);
   
  fclose(source);
  fclose(target);
}

void store_state(struct in_cmd* cmd)
{/* This functions saves the information to restart a script later */
  char tmperrn [100], tmperrn2 [100];
  char tmp_form[10];
  FILE *fptr;
  
  set_command_par_value("beam",cmd->clone, 1);

  //Change to hexdecimal format for output
  strcpy(tmp_form, float_format);
  strcpy(float_format, "A"); 
 
  char* fname     = command_par_string("file", cmd->clone);
  char* dir_name  = command_par_string("folder", cmd->clone);
  #if defined(_WIN32)
    _mkdir(dir_name);
  #elif defined(__linux__)
    mkdir(dir_name, 0700);
  #endif
  
  //The file to re-run the script
  strcpy(tmperrn2, fname);
  strcat(tmperrn2, ".madx");
  fptr = fopen(tmperrn2, "w");

  // write out the seq
  strcpy(tmperrn, fname);
  strcat(tmperrn, "_seq");
  fprintf(fptr, "call, file = %s ;\n ", tmperrn);
  fprintf(fptr,"use, sequence = %s ; \n" , current_sequ->name);
  
  // saves the macros
  strcpy(tmperrn, fname);
  strcat(tmperrn, "_macro");
  save_macros2file(tmperrn);
  fprintf(fptr, "call, file = %s ; \n",tmperrn);
  move_files(tmperrn, "", dir_name);
  
  //Save all errors
  set_selected_errors();
  if(error_esave(cmd)==1){
    strcpy(tmperrn, fname);
    strcat(tmperrn, "_errorsall");
    move_files(fname, "_errorsall", dir_name);
    fprintf(fptr, "Readmytable, file=%s, table=allerrors; \n", tmperrn);
    fprintf(fptr, "Seterr, table=%s ;\n" , "allerrors");

  }
  
  //saves the sequences
  exec_save(cmd);
  move_files(fname, "_seq", dir_name);
  
  // save and move the file created
  fclose(fptr);
  move_files(tmperrn2, "", dir_name);
  
  //Restore the old format
  strcpy(float_format, tmp_form);

  //makes a copy of the input script
  copy_input_script(strcat(dir_name, "/input_copied.madx"));
}
