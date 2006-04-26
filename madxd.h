/* should work unchanged on _win32 using Lahey */
#define advance_node          advance_node_
#define advance_to_pos        advance_to_pos_
#define augment_count         augment_count_
#define augmentcountonly      augmentcountonly_
#define char_from_table       char_from_table_  /* OB 2.4.2002 */
#define comment_to_table      comment_to_table_
#define comm_para             comm_para_
#define double_from_table     double_from_table_
#define string_from_table     string_from_table_ /* ETDA 8 nov 2004 */
#define double_to_table_row   double_to_table_row_ /* ETDA 11 nov 2004 */
#define double_table          double_table_    /* ETDA 25 aug 2004 */
#define double_to_table       double_to_table_

/* added by E. T. d'Amico on jan. 21st, 2004 */
#define interp_node           interp_node_
#define reset_interpolation   reset_interpolation_
#define embedded_twiss        embedded_twiss_
/* added by E. T. d'Amico on jan. 21stmay 19th, 2004 */
#define embedded_plot         embedded_plot_
/* end additions */

#define element_name          element_name_
#define frndm                 frndm_
#define madx                  madx_
#define madx_init             madx_init_
#define f_ctof                f_ctof_
#define get_disp0             get_disp0_
#define get_node_vector       get_node_vector_
#define get_option            get_option_
#define get_string            get_string_
#define get_title             get_title_
#define get_variable          get_variable_
#define get_vector            get_vector_
#define get_value             get_value_
#define get_version           get_version_
#define grndm                 grndm_
#define intrac                intrac_
#define mtcond                mtcond_
#define next_constraint       next_constraint_
#define next_global           next_global_
#define getnumberoftracks     getnumberoftracks_
#define next_start            next_start_
#define next_vary             next_vary_
/* RDM 20.1.2006 BEGIN jacobian strategy (match) */
#define constraint_name       constraint_name_
#define vary_name             vary_name_
/* RDM 20.1.2006 END jacobian strategy (match) */
#define node_al_errors        node_al_errors_
#define node_fd_errors        node_fd_errors_
#define node_string           node_string_
#define node_value            node_value_
#define plot_option           plot_option_
#define reset_count           reset_count_
#define restart_sequ          restart_sequ_
#define retreat_node          retreat_node_
#define sector_out            sector_out_
#define sequence_name         sequence_name_
#define set_option            set_option_
#define set_value             set_value_
#define set_variable          set_variable_
#define spec_node_value       spec_node_value_
#define store_node_value      store_node_value_
#define store_node_vector     store_node_vector_
#define string_to_table       string_to_table_
#define table_length          table_length_
#define table_org             table_org_
#define table_range           table_range_
#define track_pteigen         track_pteigen_
#define vector_to_table       vector_to_table_
#define vdot                  vdot_
#define vmod                  vmod_
#define w_ptc_create_universe   w_ptc_create_universe_
#define w_ptc_create_layout     w_ptc_create_layout_
#define w_ptc_move_to_layout    w_ptc_move_to_layout_
#define w_ptc_input             w_ptc_input_
#define w_ptc_align             w_ptc_align_
#define w_ptc_twiss             w_ptc_twiss_
#define w_ptc_normal            w_ptc_normal_
#define w_ptc_track             w_ptc_track_
#define w_ptc_start             w_ptc_start_
#define w_ptc_select            w_ptc_select_
#define w_ptc_script            w_ptc_script_
#define w_ptc_addpush           w_ptc_addpush_
#define w_ptc_end               w_ptc_end_
#define w_ptc_dumpmaps          w_ptc_dumpmaps_
#define w_ptc_trackline         w_ptc_trackline_
#define w_ptc_twiss_linac       w_ptc_twiss_linac_
#define w_ptc_setdebuglevel     w_ptc_setdebuglevel_
#define w_ptc_setaccel_method   w_ptc_setaccel_method_
#define w_ptc_setexactmis       w_ptc_setexactmis_
#define w_ptc_setradiation      w_ptc_setradiation_
#define w_ptc_setfringe         w_ptc_setfringe_
#define w_ptc_settotalpath      w_ptc_settotalpath_
#define w_ptc_settime           w_ptc_settime_
#define w_ptc_setnocavity       w_ptc_setnocavity_

#define stolower                stolower_
#define cf77flush               cf77flush_
#define select_ptc_idx          select_ptc_idx_  /* ETDA 10 nov 2004 */
#define min_order               min_order_       /* ETDA 17 nov 2004 */
#define result_from_normal      result_from_normal_ /* ETDA 11 nov 2004 */
#define make_map_table          make_map_table_ /* KZ 28.06.2005 table for maps */
#define minimum_acceptable_order minimum_acceptable_order_ /* ETDA 17 nov 2004 */

#define augmentfwarn            augmentfwarn_
/* short utility routines */
int is_operand(char c) { return (isalnum(c) || c == '_' || c == '.');}
int is_operator(char c) {return (strchr("-+*/^", c) ? 1 : 0);}
int is_expr_start(char c) {return (strchr("-+(",c) || is_operand(c));}
int mymax(int a, int b) {return (a > b ? a : b);}
int mymin(int a, int b) {return (a < b ? a : b);}
int str_pos(const char s[], char c)
{unsigned int i; for (i = 0; i < strlen(s); i++) if (s[i] == c) return i; return -1;}

/* Fortran routines called from C */
extern void dynap_(double*, double*, int*, int*, double*, double*, double*,
                   double*, double*);
/* drop first int* passed variable extern void mtgetc_(int*, double*, double*); 05.02.2005 */
extern void mtgetc_(double*, double*); /* mtgeti->mtgetc JMJ, 8/4/2003 */
extern void collect_(int*, double*, double*); /* OB 13.2.2002 */
extern void emit_(double*, double*, double*, double*, double*, double*,
                  double*, double*, double*, double*, double*, double*, double*, double*);
extern void fortinit_();
extern void getclor_(double*, double*, double*, int*);
extern void gxterm_();
extern void haveit_(double *,double *,double *,double *,int *,int *,
                    int *,double *,double *,double *,double *,double *,double *);
extern void testit_(double *,double *,double *,double *,int *,int *,
                    int *,double *,double *,double *,double *,double *,double *);
extern void svdcorr_m_(double *,double *, double *,double *,double *, double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *,int *, int *, int *);
extern void svdcorr_c_(double *,double *, double *,double *,double *, double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *,int *, int *, int *);
extern void svddec_m_(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *, int *, int *, int *);
extern void svddec_c_(double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *, int *, int *, int *);

extern void ibs_();
extern void touschek_();
extern void micit_(double *,char *,double *,double *,double *,int *,float *,
                   int *,int *,int *,int *,float *,float *,float *,float *,float *,
                   float *,float *,float *,float *,float *,int *);
extern void mtlmdf_(int*, int*, double*, int*, int*, double*, double*,
                    double*, double*, double*, double*, double*, double*,
                    double*, double*, double*, double*, double*);
extern void mtjac_(int*, int*,
                  int*, double*, double*, double*,
                   int*, int*, int*,
                   double*, int*, int*, double*, double*, double*,
                   double*, double*, double*,
                   double*, double*);
extern void mtmigr_(int*, int*, int*, double*, int*, int*, double*, double*,
                    double*, double*, double*, double*, double*, double*,
                    double*, double*, double*);
extern void mtsimp_(int*, int*, double*, int*, int*, double*, double*,
                    double*, double*, double*, double*);
extern void myindex(int*, int*, int*, int*, int*, int*, int *, int*);
extern void res_index_(int*, int*, int*, int*, int indexa[4][1000], int*);
extern void pefill_(int*);
extern void pemima_();
extern void pesopt_(int*);
extern void plotit_(int*);
extern void setup_(double *respx,double *dmat,int *im,
                   int *ic, int *nm, int*nc);
extern void soddin_(int*);
extern void survey_();
extern void tmrefe_(double*);
extern void tmrefo_(int*,double*,double*,double*);
extern void trrun_(int*,int*,double*,double*,int*,int*,
                   double*,double*,double*,double*,double*,double*,
                   double*,int*, int*, double*);
extern void twiss_(double*, double*, int*);

/* C routines called from Fortran and C */
int advance_node();
int advance_to_pos(char*, int*);
char* alias(char*);
int aperture_count(struct sequence*);
void augment_count(char*);
int char_from_table(char*, char*, int*, char*); /* OB 2.4.2002 */
void comment_to_table(char*, char*, int*);
void comm_para(char*, int*, int*, int*, int*, double*, char*, int*);
int double_from_table(char*, char*, int*, double*);
int string_from_table(char*, char*, int*, char*);
void double_to_table(char*, char*, double*);
void double_to_table_row(char*, char*, int*, double*); /* ETDA 11 nov 2004 */
int result_from_normal(char*, int*, double*); /* ETDA 11 nov 2004 */
void make_map_table(int*); /* KZ 28.06.2005 table for maps */
void element_name(char*, int*);
double frndm();
double get_aperture(struct node*, char*);
void get_disp0(double*);
void get_node_vector(char*, int*, double*);
int get_option(char*);
int get_string(char*, char*, char*);
void get_title(char*, int*);
double get_value(char*, char*);
double get_variable(char*);
int get_vector(char*, char*, double*);
void get_version(char*, int*);
double grndm();
int intrac();
int next_constraint(char*, int*, int*, double*, double*, double*, double*);
int next_global(char*, int*, int*, double*, double*, double*, double*);
int next_start(double*,double*,double*,double*,double*,double*,double*,
               double*,double*,double*,double*,double*);
/* RDM 20.1.2006 BEGIN next_vary chage definition, new func defs*/
int next_vary(char*, int*, double*, double*, double*, int*, double*);
int vary_name(char*, int*, int*);
int constraint_name(char*, int*, int*);
/* RDM 20.1.2006 END */
int node_al_errors(double*);
int node_fd_errors(double*);
void node_string(char*, char*, int*);
double node_value(char*);
void store_node_value(char* par, double* value);
double plot_option(char*);
void reset_count(char*);
int restart_sequ();
int retreat_node();
void sequence_name(char*, int*);
void set_value(char*, char*, double*);
void set_variable(char*, double*);
double spec_node_value(char*, int*);
void store_node_vector(char*, int*, double*);
void string_to_table(char*, char*, char*);
int table_length(char*);
int table_org(char*);
void table_range(char*, char*, int*);
void vector_to_table(char*, char*, int*, double*);

/* added by E. T. d'Amico */
int interp_node(int *nint);
int reset_interpolation(int *nint);
int embedded_twiss();
int select_ptc_idx(); /* 10 nov 2004 */
int minimum_acceptable_order(); /* 17 nov 2004 */
/* end additions */

/* C routines called from C */
double act_value(int, struct name_list*);
int act_special(int, char*);
int add_drifts(struct node*, struct node*);
void add_table_vars(struct name_list*, struct command_list*);
void add_to_command_list(char*, struct command*, struct command_list*, int);
void add_to_command_list_list(char*, struct command_list*,
                              struct command_list_list*);
void add_to_constraint_list(struct constraint*, struct constraint_list*);
void add_to_el_list(struct element**, int, struct el_list*, int);
void add_to_macro_list(struct macro*, struct macro_list*);
int add_to_name_list(char*, int, struct name_list*);
void add_to_node_list(struct node*, int, struct node_list*);
void add_to_sequ_list(struct sequence*, struct sequence_list*);
void add_to_table_list(struct table*, struct table_list*);
void add_to_var_list(struct variable*, struct var_list*, int);
void add_vars_to_table(struct table*);
void adjust_beam();
void all_node_pos(struct sequence*);
int attach_beam(struct sequence*);
int belongs_to_class(struct element*, char*);
char* buffer(char*);
struct in_cmd* buffered_cmd(struct in_cmd*);
void buffer_in_cmd(struct in_cmd*);
int char_cnt(char, char*);
int char_p_pos(char*, struct char_p_array*);
void check_table(char*);
struct char_p_array* clone_char_p_array(struct char_p_array*);
struct command* clone_command(struct command*);
struct command_parameter* clone_command_parameter(struct command_parameter*);
struct double_array* clone_double_array(struct double_array*);
struct element* clone_element(struct element*);
struct expression* clone_expression(struct expression*);
struct expr_list* clone_expr_list(struct expr_list*);
struct int_array* clone_int_array(struct int_array*);
struct macro* clone_macro(struct macro*);
struct name_list* clone_name_list(struct name_list*);
struct var_list* clone_var_list(struct var_list*);
struct node* clone_node(struct node*, int);
void copy_double(double*, double*, int);
void copy_name_list(struct name_list*, struct name_list*);
int cmd_match(int, char**, int*, int*);
void complete_twiss_table(struct table*);
char* compound(char*, int);
struct expression* compound_expr(struct expression*, double, char*,
                                 struct expression*, double);
void control(struct in_cmd*);
void conv_char(char*, struct int_array*);
void conv_sixtrack(struct in_cmd*);
void correct_correct(struct in_cmd*);
void correct_correct1(struct in_cmd*);
void correct_correct2(struct in_cmd*);
void correct_getorbit(struct in_cmd*);
void correct_putorbit(struct in_cmd*);
void correct_usekick(struct in_cmd*);
void correct_usemonitor(struct in_cmd*);
void correct_option(struct in_cmd* cmd);
int pro_correct_filter(int iplane, double sigcut);
int pro_correct_getcorrs(struct in_cmd* cmd);
void correct_readcorr(struct in_cmd* cmd);
void correct_setcorr(struct in_cmd* cmd);
void deco_init();
int decode_command();
int decode_par(struct in_cmd*, int, int, int, int);
struct char_array* delete_char_array(struct char_array*);
struct char_p_array* delete_char_p_array(struct char_p_array*, int);
struct command* delete_command(struct command*);
struct command_list* delete_command_list(struct command_list*);
struct command_parameter* delete_command_parameter(struct command_parameter*);
struct command_parameter_list*
delete_command_parameter_list(struct command_parameter_list*);
struct double_array* command_par_array(char*, struct command*);
struct expression* command_par_expr(char*, struct command*);
char* command_par_string(char*, struct command*);
double command_par_value(char*, struct command*);
int command_par_value2(char* parameter, struct command* cmd, double* val);
int command_par_vector(char*, struct command*, double*);
struct constraint* delete_constraint(struct constraint*);
struct constraint_list* delete_constraint_list(struct constraint_list*);
struct element* delete_element(struct element*);
struct el_list* delete_el_list(struct el_list*);
struct expression* delete_expression(struct expression*);
struct expr_list* delete_expr_list(struct expr_list*);
struct double_array* delete_double_array(struct double_array*);
struct in_cmd* delete_in_cmd(struct in_cmd*);
struct int_array* delete_int_array(struct int_array*);
struct macro* delete_macro(struct macro*);
struct name_list* delete_name_list(struct name_list*);
struct node* delete_node(struct node*);
struct node* delete_node_ring(struct node*);
struct node_list* delete_node_list(struct node_list*);
struct sequence* delete_sequence(struct sequence*);
struct sequence_list* delete_sequence_list(struct sequence_list*);
struct variable* delete_variable(struct variable*);
struct var_list* delete_var_list(struct var_list*);
struct vector_list* delete_vector_list(struct vector_list*);
struct table* delete_table(struct table*);
void disable_line(char*, struct macro_list*);
double double_from_expr(char**, int, int);
int down_unit(char*);
void dump_constraint_list(struct constraint_list*);
void dump_char_array(struct char_array*);
void dump_char_p_array(struct char_p_array*);
void dump_command(struct command*);
void dump_command_parameter(struct command_parameter*);
void dump_constraint(struct constraint*);
void dump_element(struct element*);
void dump_element_array(struct element**);
void dump_el_list(struct el_list*);
void dump_expression(struct expression*);
void dump_exp_sequ(struct sequence*, int);
void dump_in_cmd(struct in_cmd*);
void dump_int_array(struct int_array*);
void dump_macro(struct macro*);
void dump_macro_list(struct macro_list*);
void dump_name_list(struct name_list*);
void dump_node(struct node*);
void dump_sequ(struct sequence*, int);
void dump_variable(struct variable*);
void dynap_tables_create(struct in_cmd*);
double element_value(struct node*, char*);
double el_par_value(char*, struct element*);
int element_vector(struct element*, char*, double*);
void enter_element(struct in_cmd*);
void enter_elm_reference(struct in_cmd*, struct element*, int);
void enter_sequ_reference(struct in_cmd*, struct sequence*);
void enter_sequence(struct in_cmd*);
void enter_variable(struct in_cmd*);
void exec_assign(struct in_cmd*);
void exec_beam(struct in_cmd*, int);
void exec_call(struct in_cmd*);
void exec_command();
void exec_create_table(struct in_cmd*);
void exec_dump(struct in_cmd*);
void exec_dumpsequ(struct in_cmd*);
void exec_fill_table(struct in_cmd*);
void exec_help(struct in_cmd*);
void exec_macro(struct in_cmd*, int);
void exec_option();
void exec_plot(struct in_cmd*);
void exec_print(struct in_cmd*);
void exec_save(struct in_cmd*);
void exec_savebeta();
void exec_show(struct in_cmd*);
void exec_sodd(struct in_cmd*);
void exec_store_coguess(struct in_cmd*);
void expand_curr_sequ(int);
void expand_line(struct char_p_array*);
struct node* expand_node(struct node*, struct sequence*, struct sequence*,
                         double);
void expand_sequence(struct sequence*, int);
void export_comm_par(struct command_parameter*, char*);
void export_element(struct element*, struct el_list*, FILE*);
void export_elem_8(struct element*, struct el_list*, FILE*);
void export_el_def(struct element*, char*);
void export_el_def_8(struct element*, char*);
void export_el_par_8(struct command_parameter*, char*);
void export_sequence(struct sequence*, FILE*);
void export_sequ_8(struct sequence*, struct command_list*, FILE*);
void export_variable(struct variable*, FILE*);
void export_var_8(struct variable*, FILE*);
double expression_value(struct expression*, int);
void fatal_error(char*, char*);
void fill_beta0(struct command*, struct node*);
void fill_constraint_list(int, struct command*, struct constraint_list*);
void fill_elem_var_list(struct element*, struct el_list*, struct var_list*);
void fill_expr_list(char**, int, int, struct expr_list*);
void fill_expr_var_list(struct el_list*,
                        struct expression*, struct var_list*);
void fill_orbit_table(struct table*, struct table*);
void fill_par_var_list(struct el_list*,
                       struct command_parameter*, struct var_list*);
void fill_sequ_var_list(struct sequence_list*, struct el_list*,
                        struct var_list*);
void fill_twiss_header(struct table*);
void fill_twiss_header_ptc(struct table*, double);
struct command* find_command(char*, struct command_list*);
struct command_list* find_command_list(char*, struct command_list_list*);
struct element* find_element(char*, struct el_list*);
struct variable* find_variable(char*, struct var_list*);
double find_value(char*, int, char**);
int force_pos(char*);
void ftoi_array(struct double_array*, struct int_array*);
void madx();
void madx_finish();
int get_token_list(char*, char**,int);
void madx_init();
void madx_start();
void get_bracket_range(char*, char, char, int*, int*);
void get_bracket_t_range(char**, char, char, int, int, int*, int*);
void get_defined_commands();
void get_defined_constants();
struct element* get_drift(double);
char* get_new_name();
double get_node_pos(struct node*, struct sequence*);
int get_node_count(struct node*);
int get_ex_range(char*, struct sequence*, struct node**);
int get_range(char*, struct sequence*, struct node**);
double get_refpos(struct sequence*);
int get_select_ex_ranges(struct sequence*,struct command_list*,
                         struct node_list*);
int get_select_ranges(struct sequence*,struct command_list*,
                      struct node_list*);
void get_select_t_ranges(struct command_list*,
                         struct command_list*, struct table*);
int get_sub_range(char*, struct sequence*, struct node**);
int get_val_num(char*, int, int);
int square_to_colon(char*);
int get_stmt(FILE*, int);
int get_table_range(char*, struct table* t, int*);
void grow_char_array(struct char_array*);
void grow_char_array_list(struct char_array_list*);
void grow_char_p_array(struct char_p_array*);
void grow_command_list(struct command_list*);
void grow_command_list_list(struct command_list_list*);
void grow_command_parameter_list(struct command_parameter_list*);
void grow_constraint_list(struct constraint_list*);
void grow_double_array(struct double_array*);
void grow_el_list(struct el_list*);
void grow_expr_list(struct expr_list*);
void grow_in_buff_list(struct in_buff_list*);
void grow_in_cmd_list(struct in_cmd_list*);
void grow_int_array(struct int_array*);
void grow_macro_list(struct macro_list*);
void grow_name_list(struct name_list*);
void grow_node_list(struct node_list*);
void grow_sequence_list(struct sequence_list*);
void grow_table(struct table*);
void grow_table_list(struct table_list*);
void grow_var_list(struct var_list*);
void grow_vector_list(struct vector_list*);
double hidden_node_pos(char*, struct sequence*);
void init55(int);
void irngen();
int inbounds(char*, int, char**);
int in_spec_list(char*);
int int_in_array(int, int, int*);
void insert_elem(struct sequence*, struct node*);
void install_one(struct element*, char*, double, struct expression*, double);
char* join(char**, int);
char* join_b(char**, int);
int join_prefix(char*, int, char**);
double line_nodes(struct char_p_array*);
void link_in_front(struct node*, struct node*);
int loc_expr(char**, int, int, int*);
int logic_expr(int, char**);
int log_val(char*, struct command*);
void main_input(int);
struct constraint* make_constraint(int, struct command_parameter*);
struct element* make_element(char*, char*, struct command*, int);
void make_elem_node(struct element*, int);
struct expression* make_expression(int, char**);
int make_line(char*);
int make_macro(char*);
void make_occ_list(struct sequence*);
struct table* make_optics_table(struct table*);
void make_sequ_from_line(char*);
void make_sequ_node(struct sequence*, int);
char* make_string_variable(char*);
struct table* make_table(char*, char*, char**, int*, int);
void makethin(struct in_cmd*);
void match_action(struct in_cmd*);
void match_cell(struct in_cmd*);
void match_constraint(struct in_cmd*);
void match_couple(struct in_cmd*);
void match_end(struct in_cmd*);
void match_fix(struct in_cmd*);
void match_gweight(struct in_cmd*);
void match_global(struct in_cmd*);
int match_input(struct command*);  /* OB 23.1.2002 */
void match_level(struct in_cmd*);
void match_match(struct in_cmd*);
void match_prepare_varypos();
void match_rmatrix(struct in_cmd*);
void match_tmatrix(struct in_cmd*);
void match_vary(struct in_cmd*);
void match_weight(struct in_cmd*);
void mtcond(int*, int*, double*, int*);
void mtjacprint(int, int, double*);
double mult_par(char*, struct element*);
void mycpy(char*, char*);
void* mycalloc(char*, size_t, size_t);
void myfree(char*, void*);
void* mymalloc(char*, size_t);
char* mystrchr(char*, char);
void mystrcpy(struct char_array*, char*);
char* mystrstr(char*, char*);
void myrepl(char*, char*, char*, char*);
int name_list_pos(char*, struct name_list*);
struct in_buff_list* new_in_buff_list(int);
struct char_array* new_char_array(int);
struct char_array_list* new_char_array_list(int);
struct char_p_array* new_char_p_array(int);
struct command* new_command(char*, int, int, char*, char*, int, int);
struct command_list* new_command_list(char*, int);
struct command_list_list* new_command_list_list(int);
struct constraint* new_constraint(int);
struct constraint_list* new_constraint_list(int);
struct in_buffer* new_in_buffer(int);
struct in_cmd* new_in_cmd(int);
struct in_cmd_list* new_in_cmd_list(int);
struct command_parameter* new_command_parameter(char*, int);
struct command_parameter_list* new_command_parameter_list(int);
struct double_array* new_double_array(int);
struct element* new_element(char*);
struct el_list* new_el_list(int);
struct expr_list* new_expr_list(int);
struct expression* new_expression(char*, struct int_array*);
struct int_array* new_int_array(int);
struct macro* new_macro(int, int, int);
struct macro_list* new_macro_list(int);
struct name_list* new_name_list(char*, int);
struct node* new_elem_node(struct element*, int);
struct node* new_node(char*);
struct node_list* new_node_list(int);
struct sequence* new_sequence(char*, int);
struct sequence_list* new_sequence_list(int);
struct node* new_sequ_node(struct sequence*, int);
struct table* new_table(char*, char*, int, struct name_list*);
struct table_list* new_table_list(int);
struct variable* new_variable(char*, double, int, int, struct expression*,
                              char*);
struct var_list* new_var_list(int);
struct vector_list* new_vector_list(int);
int next_char(char, char**, int, int);
char next_non_blank(char*);
int next_non_blank_pos(char*);
char* noquote(char*);
void out_table(char*, struct table*, char*);
int par_present(char*, struct command*, struct command_list*);
int par_out_flag(char*, char*);
int pass_select(char*, struct command*);
int pass_select_list(char*, struct command_list*);
char* permbuff(char*);
int polish_expr(int, char**);
double polish_value(struct int_array*);
int predef_const(struct variable*);
void prepare_table_file(struct table*, struct command_list*);
void pre_split(char*, struct char_array*, int);
void print_command(struct command*);
void print_command_parameter(struct command_parameter*);
void print_global(double);
void print_rfc();
void print_table(struct table*);
void print_value(struct in_cmd*);
void pro_aperture(struct in_cmd*);
void pro_match(struct in_cmd*);
void pro_node(int, double);
void process();
void pro_correct(struct in_cmd*);
void pro_emit(struct in_cmd*);
void pro_error(struct in_cmd*);
void pro_ibs(struct in_cmd*);
void pro_touschek(struct in_cmd*);
void pro_input(char*);
void pro_sxf(struct in_cmd*);
void pro_survey(struct in_cmd*);
void pro_track(struct in_cmd*);
void pro_twiss();
void pro_ptc_twiss();
void pro_ptc_track(struct in_cmd*);
void put_info(char*, char*);
struct table* read_table(struct in_cmd*);
struct table* read_my_table(struct in_cmd*);
struct table* read_his_table(struct in_cmd*);
int remove_colon(char**, int, int);
void remove_from_command_list(char*, struct command_list*);
void remove_from_macro_list(struct macro*, struct macro_list*);
int remove_from_name_list(char*, struct name_list*);
void remove_from_node_list(struct node*, struct node_list*);
void remove_from_sequ_list(struct sequence*, struct sequence_list*);
int remove_one(struct node*);
void remove_range(char*, char*, char*);
void remove_upto(char*, char*);
void replace(char*, char, char);
void replace_lines(struct macro*, int, char**);
void replace_one(struct node*, struct element*);
void resequence_nodes(struct sequence*);
void reset_errors(struct sequence*);
void reset_sector(struct sequence*, int);
double rfc_slope();
void save_beam(struct sequence*, FILE*);
int scan_expr(int, char**);
void scan_in_cmd(struct in_cmd*);
void sector_out(double*, double*, double*, double*);
void seq_cycle(struct in_cmd*);
void seq_edit(struct in_cmd*);
void seq_edit_ex(struct sequence*);
void seq_edit_main(struct in_cmd*);
void seq_end(struct in_cmd*);
void seq_end_ex();
void seq_flatten(struct sequence*);
void seq_install(struct in_cmd*);
void seq_move(struct in_cmd*);
void seq_reflect(struct in_cmd*);
void seq_replace(struct in_cmd*);
void seq_remove(struct in_cmd*);
void set_command_par_value(char*, struct command*, double);
void set_defaults(char*);
int set_enable(char*, struct in_cmd*);
void set_new_position(struct sequence*);
void set_node_bv(struct sequence*);
void set_option(char*, int*);
void set_range(char*, struct sequence*);
void set_selected_columns(struct table*, struct command_list*);
void set_selected_elements();
void set_selected_errors();
void set_selected_rows(struct table*, struct command_list*,
                       struct command_list*);
void set_twiss_deltas(struct command*);
void set_sub_variable(char*, char*, struct in_cmd*);
void set_sector();
void show_beam(char*);
double simple_double(char**, int, int);
int simple_logic_expr(int, char**);
char* spec_join(char**, int); /* puts table() argument back for output */
int mysplit(char*, struct char_p_array*);
char* stolower(char*);  /* string to lower case in place */
void stolower_nq(char*);  /* string to lower case in place except quotes */
char* stoupper(char*);  /* string to upper case in place */
void store_beta0(struct in_cmd*);
void store_command_def(char*);
struct command_parameter* store_comm_par_def(char**, int, int);
void store_comm_par_value(char*, double, struct command*);
void store_comm_par_vector(char*, double*, struct command*);
void store_deselect(struct in_cmd*);
void store_orbit(struct command*, double*);
void store_savebeta(struct in_cmd*);
void store_select(struct in_cmd*);
void store_set(struct command*, int);
void store_threader(struct in_cmd*);
int string_cnt(char, int, char**);
char* strip(char*);
void supp_char(char, char*);
int supp_lt(char*, int);
void supp_mul_char(char, char*);
char* supp_tb(char*);
double table_value();
int table_row(struct table*, char*);
int tab_name_code(char*, char*);
void termination_handler(int);
void time_stamp(char*);
double tgrndm(double);
char* tmpbuff(char*);
void track_dynap(struct in_cmd*);
void track_end(struct in_cmd*);
void ptc_track_end();
void track_observe(struct in_cmd*);
void ptc_track_observe(struct in_cmd*);
void track_pteigen(double*);
void track_run(struct in_cmd*);
void track_ripple(struct in_cmd*);
void track_start(struct command*);
void track_tables_create(struct in_cmd*);
void track_tables_dump();
void track_track(struct in_cmd*);
void w_ptc_create_universe();
void w_ptc_create_layout();
void w_ptc_move_to_layout();
void w_ptc_input();
void w_ptc_align();
void w_ptc_twiss();
void w_ptc_normal();
void w_ptc_track();
void w_ptc_start();
void w_ptc_end();
void w_ptc_dumpmaps();
void w_ptc_trackline(int* nobspoints);
void w_ptc_twiss_linac(int* tabname);
void w_ptc_setdebuglevel(int* level);
void w_ptc_setaccel_method(int* method);
void w_ptc_setexactmis(int* boolflag);
void w_ptc_setradiation(int* boolflag);
void w_ptc_setfringe(int* boolflag);
void w_ptc_settotalpath(int* boolflag);
void w_ptc_settime(int* boolflag);
void w_ptc_setnocavity(int* boolflag);
void w_ptc_script(int* scriptname);
void w_ptc_addpush(int* tabname, int* colname, int* polinomial, int* monomial);

int twiss_input(struct command*);
void update_beam();
void update_element(struct element*, struct command*);
void update_node_constraints(struct node*, struct constraint_list*);
void update_sequ_constraints(struct sequence*, struct constraint_list*);
void update_vector(struct expr_list*, struct double_array*);
void use_sequ(struct in_cmd*);
double variable_value(struct variable*);
double vdot(int*, double*, double*);
int version_header(char*);
int v_length(char*);
char* v_format(char*);
double vmod(int*, double*);
void error(char* t1, register char* fmt, ...);
void warning(char* t1, register char* fmt, ...);
void warningOld(char*, char*);
void augmentfwarn() ;
void write_elems(struct el_list*, struct command_list*, FILE*);
void write_elems_8(struct el_list*, struct command_list*, FILE*);
void write_nice(char*, FILE*);
void write_nice_8(char*, FILE*);
void write_sequs(struct sequence_list*, struct command_list*, FILE*);
void write_table(struct table*, char*);
void write_vars(struct var_list*, struct command_list*, FILE*);
void write_vars_8(struct var_list*,struct command_list*,  FILE*);
void zero_double(double*, int);
int zero_string(char*);

/* define orbit correction routines */
void pro_correct(struct in_cmd* cmd);
int  pro_correct_getcommands(struct in_cmd* cmd);
int  pro_correct_gettables(int ip, struct in_cmd* cmd);
int  pro_correct_getorbit(struct in_cmd* cmd);
int  pro_correct_getorbit_ext(struct in_cmd* cmd);
int  pro_correct_getactive(int ip, int *nm, int *nx, int *nc, double *corvec, double *monvec,char *conm);
int  pro_correct2_gettables(int ip, struct in_cmd* cmd);
int  pro_correct2_getorbit(struct in_cmd* cmd);
int  pro_correct2_getcorrs(struct in_cmd* cmd);
int  pro_correct2_getactive(int ip, int *nm, int *nx, int *nc, double *corvec, double *monvec,char *conm);
void pro_correct_prtwiss();
void pro_correct_write_cocu_table();
void pro_correct_fill_corr_table(int ip , char *name, double old, double new);
void pro_correct2_fill_corr_table(int b, int ip , char *name, double old, double new);
void pro_correct_make_corr_table();
void pro_correct2_make_corr_table();
void pro_correct_write_results(double *monvec, double *resvec, double *corvec, int *nx, int *nc, int *nm, int imon, int icor, int ip);
void pro_correct2_write_results(double *monvec, double *resvec, double *corvec, int *nx, int *nc, int *nm, int imon, int icor, int ip);
void pro_correct_make_mon_table();
void pro_correct2_make_mon_table();
void pro_correct_fill_mon_table(int ip ,char *name, double old, double new);
void pro_correct2_fill_mon_table(int ip ,char *name, double old, double new);
double crms(double *r, int m);
float  fextim();

double* pro_correct_response_ring(int ip, int nc, int nm);
double* pro_correct2_response_ring(int ip, int nc, int nm);
double* pro_correct_response_line(int ip, int nc, int nm);

/* define utilities for orbit and error routines */
int str_from_table(char* table, char* name, int* row, char* val);
int str_from_tablet(struct table *t, char* name, int* row, char* val);

/* C wrapper to allocate memory for Fortran77 */
int c_micit(double *,char *,double *,double *,double *,int *,float,int,int,int);
void c_haveit(double *,double *,double *,double *,int *,int,int);
int  c_svddec(double *,int,int,int *);
int  c_svdcorr(double *, double *, double *, double *, int *, int , int );

/* define error routines */
double cprp(double *r, int m);
double copk(double *r, int m);
double fact(int );
void pro_error(struct in_cmd* cmd);
void error_ealign(struct in_cmd* cmd);
void error_efield(struct in_cmd* cmd);
void error_efcomp(struct in_cmd* cmd);
void error_eoption(struct in_cmd* cmd);
void error_eprint(struct in_cmd* cmd);
void error_esave(struct in_cmd* cmd);
void error_seterr(struct in_cmd* cmd);
void f_ctof(int *j, char *string, int *nel);
void pro_error_make_efield_table();

/* regular expression match routines */

struct reg_token* add_tok(char, struct reg_token*);
int char_count(char, char*);
struct reg_token* convert_pattern(char*, int, int*);
void dump_tokens(struct reg_token*);
void edit_tokens(struct reg_token*, char*, char*, int);
void fill_list(char, char, struct r_char_array*);
struct reg_token* flag(struct reg_token*, int*);
void grow_r_char_array(struct r_char_array*);
int list_count(struct r_char_array*, int, char*);
struct reg_token* make_dot(struct reg_token*);
struct reg_token* make_list(struct reg_token*, char*, int, int);
int match_all(struct reg_token*, char*);
char* match_token(struct reg_token*, char*);
int myregex(char*, char*);
void myregend(char*, struct reg_token*);
int new_comb(struct reg_token*);

/* Aperture module routines */
void aper_adj_quad(double, double, double, double*, double*);
void aper_adj_halo_si(double, double, double, double, double, double*,
                      double*, int, double*, double*);
int aper_bs(char*, double*, double*, double*, double*, int*, double*, double*);
double aper_calc(double, double, double*, double*, double*,
                 int, double*, double*, double*, double*,
                 double*, double*, int, double);
int aper_chk_inside(double, double, double*, double*, double, int);
int aper_e_d_read(char*, struct aper_e_d*, int*, char*);
int aper_external_file(char*, double*, double*);
void aper_fill_quads(double*, double*, int, int*);
void aper_header(struct table*, struct aper_node*);
void aper_intersect(double, double, double, double, double, double,
                    double, double, int, int,double*, double*);
int aper_linepar(double, double, double, double, double*, double*);
double aper_online(double, double, double, double, double, double, double);
void aper_race(double, double, double, double, double*, double*);
void aper_read_twiss(char*, int*, double*, double*, double*,
                     double*, double*, double*, double*);
int aper_rectellipse(double*, double*, double*, double*, int*, double*, double*);
void aper_surv(double*, int);
int aper_tab_search(int, struct aper_e_d*, char*, int*);
void aper_trim_ws(char*, int);
void aper_write_table(char*, double*, double*, double*, double*, double*, double*,
                      char*, double*, double*, double*, double*,
                      double*, double*, double*, double*, double*,
                      double*, double*, double*, double*, double*, char*);
struct aper_node* aperture(char*, struct node**, struct table*, int*);

/* SXF module routines */
int  all_blank(char*);
char* bpad(char*, int);
void fill_dump(FILE*, int, char*, double*, int, int);
void pro_elem_sxf(FILE*);
void put_line(FILE*, char*);
void accu_line(FILE*, char*);
void get_sxf_names();
int kl_trans(char*, char*, double*, int*);
void r_indent();
void s_indent(int);
int sxf_align_fill(int, int, int, char**, double*);
void sxf_body_fill(struct command*, int, int, int, char**, double);
int sxf_decin(char*, int);
int sxf_field_fill(int, int, int, char**, double*);
void sxf_fill_command(struct command*, int, char**);
void sxf_init();
void sxf_out();
void sxf_rtag();
void sxf_read(struct command*);
void sxf_write(struct command*, FILE*);
char* tag_spec(char*);
void reset_line(FILE*);
void write_body(FILE*);
void write_align(FILE*, struct double_array*);
void write_elend(FILE*);
void write_field(FILE*, struct double_array*);
void write_elstart(FILE*);
void cf77flush();

/*Debug level */
int debuglevel = 1;

/* Global structure variables by type (alphabetic) */
struct char_array* aux_buff;       /* temporary buffer for many purposes */
struct char_array* c_dum;
struct char_array* c_join;
struct char_array* work;
struct char_array* l_wrk;

struct char_array_list* char_buff; /* buffer for all sorts of strings */

struct char_p_array* tmp_p_array;  /* temporary buffer for splits */
struct char_p_array* tmp_l_array;  /* temporary buffer for special commands */
struct char_p_array* line_buffer;  /* buffer for line expansion */

struct command* current_beam = NULL;    /* current reference beam */
struct command* probe_beam = NULL;      /* current beam */
struct command* options = NULL;         /* current options */
struct command* plot_options = NULL;    /* current plot options */
struct command* current_error = NULL;   /* current error command */
struct command* current_correct = NULL; /* current correct command */
struct command* current_ibs = NULL;     /* current ibs command */
struct command* current_touschek = NULL;/* current touschek command */
struct command* current_survey = NULL;  /* current survey command */
struct command* current_ptc = NULL;     /* current ptc command */
struct command* current_twiss = NULL;   /* current twiss command */
struct command* current_command = NULL; /* current command clone */
struct command* current_gweight = NULL; /* current gweight clone */
struct command* current_weight = NULL;  /* current weight clone */
struct command* current_match = NULL;   /* OB 23.1.2002: current match comm. */
struct command* current_eopt  = NULL;   /* to keep eoption command */
struct command* threader_par  = NULL;   /* threader parameters */

struct command_list* beam_list;         /* list of all beam commands */
struct command_list* beta0_list;        /* list of user defined beta0s */
struct command_list* defined_commands;  /* from dictionary */
struct command_list* error_select;      /* current error select commands */
struct command_list* optics_select;     /* current optics select commands */
struct command_list* optics_list;       /* list of optics command/sequence */
struct command_list* savebeta_list;
struct command_list* seqedit_select;    /* current seqedit select commands */
struct command_list* save_select;       /* current save select commands */
struct command_list* slice_select;      /* current slice select commands */
struct command_list* stored_commands;   /* list of stored commands */
struct command_list* stored_match_var;  /* list of match vary commands */
struct command_list* stored_track_start;/* list of track start commands */
struct command_list* sector_select;     /* current sectormap select commands */

struct command_list_list* table_deselect; /* list of table deselect lists */
struct command_list_list* table_select; /* list of all table select lists */

struct constraint_list* comm_constraints; /* for each constraint command */
struct double_array* cat_doubles;    /* Polish: constant values */
struct double_array* doubles;        /* doubles buffer */
struct double_array* twiss_deltas;   /* for deltap loop in twiss command */
struct double_array* vary_vect;      /* for matching */
struct double_array* vary_dvect;     /* for matching */
struct double_array* fun_vect;       /* for matching */
struct double_array* match_work[MATCH_WORK];  /* work space for matching */

struct el_list* drift_list;
struct el_list* element_list;
struct el_list* base_type_list;
struct el_list* selected_elements;

struct in_buff_list* in;           /* list of all active input buffers */
struct in_buff_list* pro;          /* list of active processing buffers */

struct int_array* deco;       /* Polish: coded expression */
struct int_array* cat;        /* Polish: catgories */
struct int_array* d_var;      /* Polish: variable references */
struct int_array* oper;       /* Polish: operator references */
struct int_array* func;       /* Polish: function references */
struct int_array* s_range;    /* starts of ranges */
struct int_array* e_range;    /* ends of ranges */
struct int_array* sd_range;   /* starts of deselect ranges */
struct int_array* ed_range;   /* ends of deselect ranges */
struct int_array* match_i_work[MATCH_WORK];  /* int work space for matching */

struct in_cmd* this_cmd;      /* contains command just read */
struct in_cmd* local_twiss[2] = {NULL, NULL}; /* OB 1.2.2002 */
struct in_cmd* embedded_twiss_cmd = NULL;/* current plot command (ETdA 30.1.2004)*/

struct in_cmd_list* buffered_cmds;

struct macro_list* line_list;
struct macro_list* macro_list;

struct name_list* expr_chunks;
struct name_list* occ_list;
struct name_list* sxf_list;

struct node* prev_node;
struct node* current_node = NULL;
struct node* debug_node = NULL;

struct node_list* selected_ranges; /* filled by some select commands */
struct node_list* sector_ranges;   /* filled by the sectormap select command */

struct sequence* current_sequ;  /* pointer to currently used sequence */
struct sequence* edit_sequ;     /* pointer to sequence being edited */

struct sequence_list* sequences;    /* pointer to sequence list */
struct sequence_list* match_sequs;  /* pointer to sequence list for match */

struct table* aperture_table;   /* current aperture table */
struct table* ibs_table;          /* current ibs table */
struct table* touschek_table;     /* current touschek table */
struct table* summ_table;         /* current twiss summary table */
struct table* twiss_table;        /* current twiss table */
struct table* twiss_table_beam1;  /* current twiss table beam1 */
struct table* twiss_table_beam2;  /* current twiss table beam2 */
struct table* map_table;          /* added for twiss_input_table */
struct table_list* table_register;/* added by kzhang 26/06/2005 */


/* E. T. d'Amico 2 feb 2004 */
struct table* embedded_twiss_table;        /* current twiss table */
/* end additions */
/* E. T. d'Amico 5 nov 2004 */
struct table* normal_results;     /* ptc table containing the selected high order functions (such as dx,qx,anhx etc.) */
/* end additions */

struct table* survey_table;       /* current survey table */
struct table* corr_table;         /* corrector table after orbit correction */
struct table* corr_table1;        /* corrector table after orbit correction, beam 1 for two rings */
struct table* corr_table2;        /* corrector table after orbit correction, beam 2 for two rings */
struct table* mon_table;          /* monitor table after orbit correction */
struct table* orbit_table;        /* orbit positions at monitors */
struct table* sodd_table_70;      /* sodd output table detune_1_end */
struct table* sodd_table_71;      /* sodd output table detune_1_all */
struct table* sodd_table_72;      /* sodd output table detune_2_end */
struct table* sodd_table_73;      /* sodd output table detune_2_all */
struct table* sodd_table_74;      /* sodd output table distort_1_F_end */
struct table* sodd_table_75;      /* sodd output table distort_1_H_end */
struct table* sodd_table_76;      /* sodd output table distort_1_F_all */
struct table* sodd_table_77;      /* sodd output table distort_1_H_all */
struct table* sodd_table_78;      /* sodd output table distort_2_F_end */
struct table* sodd_table_79;      /* sodd output table distort_2_F_all */
struct table* target_table = NULL;       /* current target table */
struct table* model_table = NULL;        /* current model table */
struct table* orbin_table = NULL;        /* current orbit table */


struct table_list* optics_tables; /* contains optics tables from last twiss */
struct table_list* table_register; /* contains all tables */

struct variable* current_variable = NULL; /* set by act_value (table access) */
struct var_list* variable_list;

struct orb_cor*  correct_orbit;   /* information and links for orbit correction */
struct orb_cor2* correct_orbit1;  /* information and links for orbit correction */
struct orb_cor2* correct_orbit2;  /* information and links for orbit correction */
struct orb_cor2* correct_orbit12; /* information and links for orbit correction */

double corrl;                  /* global limit for orbit corrector strength  */

struct table* efield_table;    /* field errors in table form  */
FILE* fddata;
FILE* fcdata;
FILE* ftdata;
FILE* fgdata;

FILE* debug_file;              /* for debug output */
FILE* stamp_file;              /* for debug output */
FILE* out_file;                /* for table output */
FILE* prt_file;                /* for echo output */
FILE* sec_file = NULL;         /* for sector output */
FILE* tab_file;                /* for table input */

/* Global simple variables by type */

char quote;                       /* current open single or double quote */
char tmp_key[NAME_L],
  int_format[20],             /* current integer format */
  float_format[20],           /* current float format */
  string_format[20];          /* current string format */
char var_form[1000];             /* buffer for the user-controlled formats */
char blank[] = "    ";
char none[] = "none";
char myversion[] = "MAD-X 3.02.29";
char code_mod_date[] = "Code Modification Date: 25.04.2006";
char one_string[] = "1";
char aptwfile[FNAME_L] = "dummy"; /* IW 02.12.2004 */
char* aux_char_pt;               /* for debug purposes */
char* exx;
char* current_link_group;
char* current_range;             /* currently used range, or NULL */
char* title = NULL;
char* match_seqs[2];             /* OB 23.1.2002   */
char* match_beta[2];             /* OB 23.1.2002   */
char* match_range[2];            /* HG 12.11.2002   */
char* myfree_caller = none;      /* for routine name in myfree calls */
char* track_filename;            /* track module file name start */
char* track_fileext;             /* track module file name extension */
/* E. T. d'Amico 11 june 2004 */
char track_plot_filename[NAME_L] = "madx_track";            /* plot module: output postscript file name in track mode */
/* end additions */

double pi, twopi, degrad, raddeg, e, clight, hbar;
double penalty;
double match_tol;
double orbit0[6];
double disp0[6];
double sxf_suml = 0;
double track_deltap=0;
double oneturnmat[36];
/* E. T. d'Amico 13 may 2004 */
double fintx_plot;              /* to save the value of fintx for the reset_interpolation routine */
/* end additions */

const double zero = 0;
const double one = 1;
const double two = 2;
const double ten_p_3 = 1.e3;
const double ten_p_6 = 1.e6;
const double ten_p_9 = 1.e9;
const double ten_p_12 = 1.e12;
const double ten_m_3 = 1.e-3;
const double ten_m_6 = 1.e-6;
const double ten_m_9 = 1.e-9;
const double ten_m_12 = 1.e-12;
const double ten_m_15 = 1.e-15;
const double ten_m_16 = 1.e-16;
const double ten_m_19 = 1.e-19;

int add_error_opt = 0;      /* ADD error option, set with eoption */
int warn_numb = 0;           /* Number of warnings */
int warn_numbf = 0;           /* Number of warnings from fortran*/
/* E. T. d'Amico 25 feb 2004 */
int rbend = 0;              /* flag (= 1 when the element is a rectangular bending magnet) */
/* E. T. d'Amico 13 may 2004 */
int embedded_flag = 0;              /* flag (= 1 when entering routine pro_embedded_twiss, 0 at exit) */
/* E. T. d'Amico 17 nov 2004 */
int min_order = 1;      /* minimum required order */
/* end additions */
int print_correct_opt = 1;  /* PRINT options for orbit correction */
int debug_correct_opt = 0;  /* DEBUG options for orbit correction */
int assign_start = 0;       /* flag for multiple assign statements */
int aux_count = 0;          /* for debug purposes */
int beam_info = -1;         /* flag to print beam information once */
int c_range_end;            /* node count of current range end */
int c_range_start;          /* node count of current range start */
int curr_obs_points;        /* current number of observation points */
int current_calls = 0;      /* call counter in match */
int current_call_lim = 0;   /* current call limit in match */
int current_const = 0;      /* current constraint number in match */
int default_beam_saved = 0; /* flag to avoid multiple save of default beam */
int edit_is_on = 0;         /* != 0 if inside current sequence edit */
int final_message = 0;      /* set to 1 when end message written */
int group_is_on = 0;        /* true when inside group */
int guess_flag = 0;         /* != 0 if coguess read */
int in_stop = 0;            /* input buffer stop flag */
int inbuf_level = 0;        /* input buffer level */
int init_warn = 1;          /* intialisation warning level */
int interactive;            /* non-zero if interactive */
int irn_rand[NR_RAND];      /* for random generator */
int keep_tw_print;          /* previous twiss print flag (match) */
int loop_cnt = 0;           /* used to detect infinite loops */
int match_calls = 0;        /* command call limit in match */
int match_is_on = 0;        /* true when inside match command */
int match_num_beta = 0;     /* OB 23.1.2002 */
int match_num_range = 0;    /* HG 12.11.2002 */
int match_num_seqs = 0;     /* OB 23.1.2002 */
int mig_strategy;           /* migrad strategy (match) */
int jac_strategy;           /* RDM 24.8.2005 jacobian strategy (match) */
int jac_repeat;             /* RDM 24.8.2005 jacobian repeat (match) */
double jac_cool;            /* RDM 24.8.2005 jacobian cool factor (match) */
double jac_balance;         /* RDM 24.8.2005 jacobian balance cool factor (match) */
double jac_random;         /* RDM 24.8.2005 jacobian random factor (match) */
int jac_bisec;             /* RDM 16.3.2006 jacobian bisec factor (match) */
int new_name_count = 0;     /* to make internal names */
int next_rand = 0;          /* for random generator */
int plots_made = 0;         /* set to 1 if plots are made */
int polish_cnt = 0;         /* used to detect infinite loops */
int print_match_summary = 0;/* OB 6.3.2002:
                               activate the print option in the
                               'mtgeti' and 'collect' routines (mtgeti->mtgetc JMJ, 8/4/2003)    */
int quote_toggle = 0;       /* for quote strings on input */
int return_flag = 0;        /* 1 when "return" read */
int scrap_count = 0;        /* running counter to make things unique */
int seqedit_install = 0;    /* counter for seqedit installs */
int seqedit_move = 0;       /* counter for seqedit moves */
int seqedit_remove = 0;     /* counter for seqedit removes */
int sequ_is_on = 0;         /* != 0 if inside current sequence decl. */
int stamp_flag = 0;         /* checks for double delete when != 0 */
int start_cnt = 0;          /* counter for start commands */
int start_var = 0;          /* start of variables after predefined constants */
int total_const = 0;        /* total no. of constraints in match */
int total_vars = 0;         /* total no. of variables in match */
int track_is_on = 0;        /* true when inside track command */
int track_start_cnt = 0;    /* counter for track start commands */
int twiss_success = 0;      /* set by twiss module to 1 if OK */
int use_count = 0;          /* incremented by 1 every time use is executed */
int vary_cnt = 0;           /* counter for vary commands */
int watch_flag = 0;         /* produces debug output when != 0 */

int           na_err,              /* current no. of alignment errors */
  nf_err,              /* current no. of field errors */
  indent = 0,          /* current indentation count */
  b_level = 0,         /* current brace level */
  sxf_elem_cnt = 0,    /* element count */
  tag_flag = 0,        /* if > 0, tag = parent name written */
  tag_cnt = 0,         /* if > 0, tag = specified type code
                          written for selected types only */
  sxf_align_cnt = 0,       /* element with align errors count */
  sxf_field_cnt = 0,       /* element with field errors count */
  stop_flag = 0,           /* 1 if stop condition */
  occnt_add = 0,       /* flag for element name modification */
  b_indent[100],       /* list of indents */
  add_indent[] = {1, 2, 2, 4, 7, 7, 7, 7, 7, 7};

double        sequ_length,         /* length of  sequence */
  sequ_start,
  sequ_end,
  guess_orbit[6],
  al_errors[ALIGN_MAX],
  fd_errors[FIELD_MAX];

char          line[LINE_MAX],
  tag_type[MAX_TAG][16],
  tag_code[MAX_TAG][16];

time_t last_time;
time_t start_time;

/*Piotr Skowronski (CERN)*/
#define gettrack gettrack_
#define deletetrackstrarpositions deletetrackstrarpositions_

double** trackstrarpositions = 0x0;/* two dimensional array with track positions*/

int  gettrack(int* n, double* x,double* px,double* y,double* py,double* t,double* pt);
int  copytrackstoarray();
void deletetrackstrarpositions();

/*Riccardo de Maria (CERN)*/
void match2_match(struct in_cmd*);
void match2_end(struct in_cmd*);
void match2_macro(struct in_cmd*);
void match2_constraint(struct in_cmd*);
char* match2_macro_name[10];
char* match2_cons_name[10][30];
double match2_cons_value[10][30];
double match2_cons_value_rhs[10][30];
double match2_cons_value_lhs[10][30];
double match2_cons_weight[10][30];
char match2_cons_sign[10][30];
int match2_cons_curr[3];
struct expression* match2_cons_rhs[10][30];
struct expression* match2_cons_lhs[10][30];

/* end of definitions */

