#include "madx.h"

/* Temporary file: global variables
   these variables will be split over their respective modules...
*/

int  debuglevel = 1;

/* Global structure variables by type (alphabetic) */

struct char_array* aux_buff;       /* temporary buffer for many purposes */
struct char_array* c_dum;
struct char_array* c_join;
struct char_array* work;
struct char_array* l_wrk;

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
struct command* current_match = NULL;   /* current match comm. */
struct command* current_eopt = NULL;    /* to keep eoption command */
struct command* threader_par = NULL;    /* threader parameters */

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
struct command_list* interp_select;     /* current interpolate select commands */

struct command_list_list* table_deselect; /* list of table deselect lists */
struct command_list_list* table_select;   /* list of all table select lists */

struct constraint_list* comm_constraints; /* for each constraint command */
struct double_array* cat_doubles;    /* Polish: constant values */
struct double_array* doubles;        /* doubles buffer */
struct double_array* twiss_deltas;   /* for deltap loop in twiss command */
struct double_array* vary_vect;      /* for matching */
struct double_array* vary_dvect;     /* for matching */
struct double_array* fun_vect;       /* for matching */
struct double_array* match_work[MATCH_WORK];/* work space for matching */

struct el_list* drift_list;
struct el_list* element_list;
struct el_list* base_type_list;
struct el_list* selected_elements;

struct in_buff_list* in;      /* list of all active input buffers */
struct in_buff_list* pro;     /* list of active processing buffers */

struct int_array* deco;       /* Polish: coded expression */
struct int_array* cat;        /* Polish: catgories */
struct int_array* d_var;      /* Polish: variable references */
struct int_array* oper;       /* Polish: operator references */
struct int_array* func;       /* Polish: function references */
struct int_array* s_range;    /* starts of ranges */
struct int_array* e_range;    /* ends of ranges */
struct int_array* sd_range;   /* starts of deselect ranges */
struct int_array* ed_range;   /* ends of deselect ranges */
struct int_array* match_i_work[MATCH_WORK];/* int work space for matching */

struct in_cmd* this_cmd;      /* contains command just read */
struct in_cmd* local_twiss[2] = {NULL, NULL};
struct in_cmd* embedded_twiss_cmd = NULL;/* current plot command */

struct in_cmd_list* buffered_cmds;

struct macro_list* line_list;
struct macro_list* macro_list;

struct name_list* expr_chunks;
struct name_list* occ_list;
struct name_list* sxf_list;

struct node* prev_node;
struct node* current_node = NULL;
struct node* debug_node = NULL;

struct node_list* selected_ranges;/* filled by some select commands */
struct node_list* sector_ranges;  /* filled by the sectormap select command */

struct sequence* current_sequ;    /* pointer to currently used sequence */
struct sequence* edit_sequ;       /* pointer to sequence being edited */

struct sequence_list* sequences;  /* pointer to sequence list */
struct sequence_list* match_sequs;/* pointer to sequence list for match */

struct table* aperture_table;     /* current aperture table */
struct table* ibs_table;          /* current ibs table */
struct table* touschek_table;     /* current touschek table */
struct table* summ_table;         /* current twiss summary table */
struct table* twiss_table;        /* current twiss table */
struct table* twiss_table_beam1;  /* current twiss table beam1 */
struct table* twiss_table_beam2;  /* current twiss table beam2 */
struct table* twiss_sector_table; /* used for sectormap */
struct table* ptc_twiss_summary_table;/* holds summary data after one turn */
struct table* map_table;          /* added for twiss_input_table */
struct table_list* table_register;
struct table_list* moments_tables = 0x0;/* tables for moments */

struct table* embedded_twiss_table;/* current twiss table */
struct table* normal_results;     /* ptc table containing the selected high order functions (such as dx,qx,anhx etc.) */

struct table* errors_dipole;
struct table* errors_field;
struct table* errors_total;
struct table* errors_read; /* table needed for IO of errors with PTC */

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
struct table* sodd_table_74;      /* sodd output table distort_1_f_end */
struct table* sodd_table_75;      /* sodd output table distort_1_h_end */
struct table* sodd_table_76;      /* sodd output table distort_1_f_all */
struct table* sodd_table_77;      /* sodd output table distort_1_h_all */
struct table* sodd_table_78;      /* sodd output table distort_2_f_end */
struct table* sodd_table_79;      /* sodd output table distort_2_h_end */
struct table* target_table = NULL;/* current target table */
struct table* model_table = NULL; /* current model table */
struct table* orbin_table = NULL; /* current orbit table */


struct table_list* optics_tables; /* contains optics tables from last twiss */
struct table_list* table_register;/* contains all tables */

struct table_list_list* all_table_lists; /* all table lists are entered here */

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

struct char_p_array* sdds_pat; /* array for selected sdds patterns */

FILE* debug_file;              /* for debug output */
FILE* stamp_file;              /* for debug output */
FILE* out_file;                /* for table output */
FILE* prt_file;                /* for echo output */
FILE* sec_file = NULL;         /* for sector output in "embedded" twiss */
FILE* tab_file;                /* for table input */

/* Global simple variables by type */

char  quote;                      /* current open single or double quote */
char  int_format[20],             /* current integer format */
      float_format[20],           /* current float format */
      string_format[20];          /* current string format */
char  blank[] = "    ";
char  none[] = "none";
char  one_string[] = "1";
// 2015-Jul-31  11:41:59  ghislain: aperture twiss file for output of twiss table ! not needed
//char  aptwfile[FNAME_L] = "dummy";
char* aux_char_pt;               /* for debug purposes */
char* exx;
char* current_link_group;
char* current_range;             /* currently used range, or NULL */
char* title = NULL;
char* match_seqs[2];
char* match_beta[2];
char* match_range[2];
char* track_filename;            /* track module file name start */
char* track_fileext;             /* track module file name extension */
char  track_plot_filename[NAME_L] = "madx_track"; /* plot module: output postscript file name in track mode */

double pi, twopi, degrad, raddeg, e, clight, hbar;
double penalty;
double match_tol;
double orbit0[6];
double disp0[6];
double sxf_suml = 0;
double track_deltap=0;
double oneturnmat[36];

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
int embedded_flag = 0;      /* flag (= 1 when entering routine pro_embedded_twiss, 0 at exit) */
int min_order = 1;          /* minimum required order */
int print_correct_opt = 1;  /* PRINT options for orbit correction */
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
int keep_tw_print;          /* previous twiss print flag (match) */
// int loop_cnt = 0;           /* used to detect infinite loops ; removed 2014-Mar-20  16:20:13  ghislain */
int match_calls = 0;        /* command call limit in match */
int match_is_on = 0;        /* true when inside match command */
int chrom_match = 0;        /* true when the match summary*/
int match_num_beta = 0;
int match_num_range = 0;
int match_num_seqs = 0;
int mig_strategy;           /* migrad strategy (match) */
int jac_strategy;           /* jacobian strategy (match) */
int jac_repeat;             /* jacobian repeat (match) */
double jac_cool;            /* jacobian cool factor (match) */
double jac_balance;         /* jacobian balance cool factor (match) */
double jac_random;          /* jacobian random factor (match) */
int jac_bisec;              /* jacobian bisec factor (match) */
double jac_cond;            /* jacobian svd cond. num (match) */
int new_name_count = 0;     /* to make internal names */
int next_rand = 0;          /* for random generator */
int plots_made = 0;         /* set to 1 if plots are made */
int polish_cnt = 0;         /* used to detect infinite loops */
int print_match_summary = 0;/* activate the print option in the
                               'mtgeti' and 'collect' routines (mtgeti->mtgetc) */
int quote_toggle = 0;       /* for quote strings on input */
int return_flag = 0;        /* 1 when "return" read */
int scrap_count = 0;        /* running counter to make things unique */
int seqedit_install = 0;    /* counter for seqedit installs */
int seqedit_move = 0;       /* counter for seqedit moves */
int seqedit_remove = 0;     /* counter for seqedit removes */
int seqedit_replace = 0;    /* counter for seqedit replaces -- 2014-Jul-01  13:28:05  ghislain */
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

int na_err,                 /* current no. of alignment errors */
    nf_err,                 /* current no. of field errors */
    indent = 0,             /* current indentation count */
    b_level = 0,            /* current brace level */
    sxf_elem_cnt = 0,       /* element count */
    tag_flag = 0,           /* if > 0, tag = parent name written */
    tag_cnt = 0,            /* if > 0, tag = specified type code
                               written for selected types only */
    sxf_align_cnt = 0,      /* element with align errors count */
    sxf_field_cnt = 0,      /* element with field errors count */
    stop_flag = 0,          /* 1 if stop condition */
    occnt_add = 0,          /* flag for element name modification */
    b_indent[100],          /* list of indents */
    add_indent[] = {1, 2, 2, 4, 7, 7, 7, 7, 7, 7};

double
    guess_orbit[6],
    al_errors[ALIGN_MAX],
    fd_errors[FIELD_MAX];

char
    line[MADX_LINE_MAX],
    tag_type[MAX_TAG][16],
    tag_code[MAX_TAG][16];

time_t last_time;
time_t start_time;

char filenames[100][500];
int  currentline[100];

double** trackstrarpositions = 0x0; /* two dimensional array with track positions*/

