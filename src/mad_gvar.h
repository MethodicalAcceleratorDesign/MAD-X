#ifndef MAD_GVAR_H
#define MAD_GVAR_H

/* Temporary file: global variables
   these variables will be split over their respective modules...
*/

extern int debuglevel;

/* Global structure variables by type (alphabetic) */

extern struct char_array* aux_buff;       /* temporary buffer for many purposes */
extern struct char_array* c_dum;
extern struct char_array* c_join;
extern struct char_array* work;
extern struct char_array* l_wrk;

extern struct char_p_array* tmp_p_array;  /* temporary buffer for splits */
extern struct char_p_array* tmp_l_array;  /* temporary buffer for special commands */
extern struct char_p_array* line_buffer;  /* buffer for line expansion */

extern struct command* current_beam;    /* current reference beam */
extern struct command* probe_beam;      /* current beam */
extern struct command* options;         /* current options */
extern struct command* plot_options;    /* current plot options */
extern struct command* current_error;   /* current error command */
extern struct command* current_correct; /* current correct command */
extern struct command* current_ibs;     /* current ibs command */
extern struct command* current_touschek;/* current touschek command */
extern struct command* current_survey;  /* current survey command */
extern struct command* current_ptc;     /* current ptc command */
extern struct command* current_twiss;   /* current twiss command */
extern struct command* current_command; /* current command clone */
extern struct command* current_gweight; /* current gweight clone */
extern struct command* current_weight;  /* current weight clone */
extern struct command* current_match;   /* current match comm. */
extern struct command* current_eopt;    /* to keep eoption command */
extern struct command* threader_par;    /* threader parameters */
extern struct command* current_format_f;   /* current formats */

extern struct command_list* beam_list;         /* list of all beam commands */
extern struct command_list* beta0_list;        /* list of user defined beta0s */
extern struct command_list* defined_commands;  /* from dictionary */
extern struct command_list* error_select;      /* current error select commands */
extern struct command_list* optics_select;     /* current optics select commands */
extern struct command_list* optics_list;       /* list of optics command/sequence */
extern struct command_list* savebeta_list;
extern struct command_list* seqedit_select;    /* current seqedit select commands */
extern struct command_list* save_select;       /* current save select commands */
extern struct command_list* slice_select;      /* current slice select commands */
extern struct command_list* stored_commands;   /* list of stored commands */
extern struct command_list* stored_match_var;  /* list of match vary commands */
extern struct command_list* stored_track_start;/* list of track start commands */
extern struct command_list* sector_select;     /* current sectormap select commands */
extern struct command_list* interp_select;     /* current interpolate select commands */

extern struct command_list_list* table_deselect; /* list of table deselect lists */
extern struct command_list_list* table_select; /* list of all table select lists */

extern struct constraint_list* comm_constraints; /* for each constraint command */
extern struct double_array* cat_doubles;    /* Polish: constant values */
extern struct double_array* doubles;        /* doubles buffer */
extern struct double_array* twiss_deltas;   /* for deltap loop in twiss command */
extern struct double_array* vary_vect;      /* for matching */
extern struct double_array* vary_dvect;     /* for matching */
extern struct double_array* fun_vect;       /* for matching */
extern struct double_array* match_work[MATCH_WORK];/* work space for matching */

extern struct el_list* element_list;        /* Explicitly defined elems. No implicit drifts! */
extern struct el_list* base_type_list;
extern struct el_list* selected_elements;

extern struct expression* backup_expr;

extern struct in_buff_list* in;      /* list of all active input buffers */
extern struct in_buff_list* pro;     /* list of active processing buffers */

extern struct int_array* deco;       /* Polish: coded expression */
extern struct int_array* cat;        /* Polish: catgories */
extern struct int_array* d_var;      /* Polish: variable references */
extern struct int_array* oper;       /* Polish: operator references */
extern struct int_array* func;       /* Polish: function references */
extern struct int_array* s_range;    /* starts of ranges */
extern struct int_array* e_range;    /* ends of ranges */
extern struct int_array* sd_range;   /* starts of deselect ranges */
extern struct int_array* ed_range;   /* ends of deselect ranges */
extern struct int_array* match_i_work[MATCH_WORK];  /* int work space for matching */

extern struct in_cmd* this_cmd;      /* contains command just read */
extern struct in_cmd* local_twiss[2];
extern struct in_cmd* embedded_twiss_cmd; /* current plot command */

extern struct in_cmd_list* buffered_cmds;

extern struct macro_list* line_list;
extern struct macro_list* macro_list;

extern struct name_list* expr_chunks;
extern struct name_list* occ_list;
extern struct name_list* sxf_list;

extern struct node* prev_node;
extern struct node* current_node;
extern struct node* debug_node;

extern struct node_list* selected_ranges; /* filled by some select commands */
extern struct node_list* sector_ranges;   /* filled by the sectormap select command */

extern struct sequence* current_sequ;  /* pointer to currently used sequence */
extern struct sequence* edit_sequ;     /* pointer to sequence being edited */

extern struct sequence_list* sequences;    /* pointer to sequence list */
extern struct sequence_list* match_sequs;  /* pointer to sequence list for match */

extern struct table* aperture_table;     /* current aperture table */
extern struct table* ibs_table;          /* current ibs table */
extern struct table* touschek_table;     /* current touschek table */
extern struct table* summ_table;         /* current twiss summary table */
extern struct table* twiss_table;        /* current twiss table */
extern struct table* twiss_table_beam1;  /* current twiss table beam1 */
extern struct table* twiss_table_beam2;  /* current twiss table beam2 */
extern struct table* twiss_sector_table; /* used for sectormap */
extern struct table* ptc_twiss_summary_table; /* holds summary data after one turn */
extern struct table* map_table;          /* added for twiss_input_table */
extern struct table_list* table_register;
extern struct table_list* moments_tables;/* tables for moments */

extern struct table* embedded_twiss_table;/* current twiss table */
extern struct table* normal_results;     /* ptc table containing the selected high order functions (such as dx,qx,anhx etc.) */

extern struct table* errors_dipole;
extern struct table* errors_field;
extern struct table* errors_total;
extern struct table* errors_read;        /* table needed for IO of errors with PTC */

extern struct table* survey_table;       /* current survey table */
extern struct table* corr_table;         /* corrector table after orbit correction */
extern struct table* corr_table1;        /* corrector table after orbit correction, beam 1 for two rings */
extern struct table* corr_table2;        /* corrector table after orbit correction, beam 2 for two rings */
extern struct table* mon_table;          /* monitor table after orbit correction */
extern struct table* orbit_table;        /* orbit positions at monitors */
extern struct table* sodd_table_70;      /* sodd output table detune_1_end */
extern struct table* sodd_table_71;      /* sodd output table detune_1_all */
extern struct table* sodd_table_72;      /* sodd output table detune_2_end */
extern struct table* sodd_table_73;      /* sodd output table detune_2_all */
extern struct table* sodd_table_74;      /* sodd output table distort_1_f_end */
extern struct table* sodd_table_75;      /* sodd output table distort_1_h_end */
extern struct table* sodd_table_76;      /* sodd output table distort_1_f_all */
extern struct table* sodd_table_77;      /* sodd output table distort_1_h_all */
extern struct table* sodd_table_78;      /* sodd output table distort_2_f_end */
extern struct table* sodd_table_79;      /* sodd output table distort_2_h_end */
extern struct table* target_table;       /* current target table */
extern struct table* model_table;        /* current model table */
extern struct table* orbin_table;        /* current orbit table */


extern struct table_list* optics_tables; /* contains optics tables from last twiss */
extern struct table_list* table_register; /* contains all tables */

extern struct table_list_list* all_table_lists; /* all table lists are entered here */

extern struct variable* current_variable; /* set by act_value (table access) */
extern struct var_list* variable_list;

extern struct orb_cor*  correct_orbit;   /* information and links for orbit correction */
extern struct orb_cor2* correct_orbit1;  /* information and links for orbit correction */
extern struct orb_cor2* correct_orbit2;  /* information and links for orbit correction */
extern struct orb_cor2* correct_orbit12; /* information and links for orbit correction */

extern double corrl;                  /* global limit for orbit corrector strength  */

extern struct table* efield_table;    /* field errors in table form  */
extern FILE* fddata;
extern FILE* fcdata;
extern FILE* ftdata;
extern FILE* fgdata;

extern struct char_p_array* sdds_pat; /* array for selected sdds patterns */

extern FILE* debug_file;              /* for debug output */
extern FILE* stamp_file;              /* for debug output */
extern FILE* out_file;                /* for table output */
extern FILE* prt_file;                /* for echo output */
extern FILE* sec_file;                /* for sector output in "embedded" twiss */
extern FILE* tab_file;                /* for table input */

/* Global simple variables by type */

extern char quote;                      /* current open single or double quote */
extern char int_format[20],             /* current integer format */
            float_format[20],           /* current float format */
            string_format[20];          /* current string format */
extern char blank[];
extern char none[];
extern char one_string[];
// 2015-Jul-31  11:41:59  ghislain: aperture twiss file for output of twiss table ! not needed
//extern char aptwfile[FNAME_L];
extern char* aux_char_pt;               /* for debug purposes */
extern char* exx;
extern char* current_link_group;
extern char* current_range;             /* currently used range, or NULL */
extern char* title;
extern char* match_seqs[2];
extern char* match_beta[2];
extern char* match_range[2];
extern char* track_filename;            /* track module file name start */
extern char* track_fileext;             /* track module file name extension */
extern char  track_plot_filename[NAME_L];/* plot module: output postscript file name in track mode */

extern double pi, twopi, degrad, raddeg, e, clight, hbar;
extern double penalty;
extern double match_tol;
extern double orbit0[6];
extern double disp0[6];
extern double sxf_suml;
extern double track_deltap;
extern double oneturnmat[36];

extern const double zero;
extern const double one;
extern const double two;
extern const double ten_p_3;
extern const double ten_p_6;
extern const double ten_p_9;
extern const double ten_p_12;
extern const double ten_m_3;
extern const double ten_m_6;
extern const double ten_m_9;
extern const double ten_m_12;
extern const double ten_m_15;
extern const double ten_m_16;
extern const double ten_m_19;

extern int add_error_opt;          /* ADD error option, set with eoption */
extern int backup_type;

extern int embedded_flag;          /* flag (= 1 when entering routine pro_embedded_twiss, 0 at exit) */
extern int min_order;              /* minimum required order */
extern int print_correct_opt;      /* PRINT options for orbit correction */
extern int assign_start;           /* flag for multiple assign statements */
extern int aux_count;              /* for debug purposes */
extern int beam_info;              /* flag to print beam information once */
extern int curr_obs_points;        /* current number of observation points */
extern int current_calls;          /* call counter in match */
extern int current_call_lim;       /* current call limit in match */
extern int current_const;          /* current constraint number in match */
extern int default_beam_saved;     /* flag to avoid multiple save of default beam */
extern int edit_is_on;             /* != 0 if inside current sequence edit */
extern int final_message;          /* set to 1 when end message written */
extern int group_is_on;            /* true when inside group */
extern int guess_flag;             /* != 0 if coguess read */
extern int in_stop;                /* input buffer stop flag */
extern int inbuf_level;            /* input buffer level */
extern int init_warn;              /* intialisation warning level */
extern int interactive;            /* non-zero if interactive */
extern int keep_tw_print;          /* previous twiss print flag (match) */
// extern int loop_cnt;               /* used to detect infinite loops ; removed 2014-Mar-20  16:20:13  ghislain */
extern int match_calls;            /* command call limit in match */
extern int match_is_on;            /* true when inside match command */
extern int chrom_match;            /* true when the match summary*/
extern int match_num_beta;
extern int match_num_range;
extern int match_num_seqs;
extern int mig_strategy;           /* migrad strategy (match) */
extern int jac_strategy;           /* jacobian strategy (match) */
extern int jac_repeat;             /* jacobian repeat (match) */
extern double jac_cool;            /* jacobian cool factor (match) */
extern double jac_balance;         /* jacobian balance cool factor (match) */
extern double jac_random;          /* jacobian random factor (match) */
extern int jac_bisec;              /* jacobian bisec factor (match) */
extern double jac_cond;            /* jacobian svd cond. num (match) */
extern int new_name_count;         /* to make internal names */
extern int next_rand;              /* for random generator */
extern int plots_made;             /* set to 1 if plots are made */
extern int polish_cnt;             /* used to detect infinite loops */
extern int print_match_summary;    /* activate the print option in the
                                      'mtgeti' and 'collect' routines (mtgeti->mtgetc) */
extern int quote_toggle;           /* for quote strings on input */
extern int return_flag;            /* 1 when "return" read */
extern int scrap_count;            /* running counter to make things unique */
extern int seqedit_install;        /* counter for seqedit installs */
extern int seqedit_move;           /* counter for seqedit moves */
extern int seqedit_remove;         /* counter for seqedit removes */
extern int seqedit_replace;        /* counter for seqedit replaces -- 2014-Jul-01  13:27:11  ghislain */
extern int sequ_is_on;             /* != 0 if inside current sequence decl. */
extern int stamp_flag;             /* checks for extern double delete when != 0 */
extern int start_cnt;              /* counter for start commands */
extern int start_var;              /* start of variables after predefined constants */
extern int total_const;            /* total no. of constraints in match */
extern int total_vars;             /* total no. of variables in match */
extern int track_is_on;            /* true when inside track command */
extern int track_start_cnt;        /* counter for track start commands */
extern int twiss_success;          /* set by twiss module to 1 if OK */
extern int use_count;              /* incremented by 1 every time use is executed */
extern int vary_cnt;               /* counter for vary commands */
extern int watch_flag;             /* produces debug output when != 0 */

extern int na_err,                 /* current no. of alignment errors */
           nf_err,                 /* current no. of field errors */
           indent,                 /* current indentation count */
           b_level,                /* current brace level */
           sxf_elem_cnt,           /* element count */
           tag_flag,               /* if > 0, tag = parent name written */
           tag_cnt,                /* if > 0, tag = specified type code
                                      written for selected types only */
           sxf_align_cnt,          /* element with align errors count */
           sxf_field_cnt,          /* element with field errors count */
           stop_flag,              /* 1 if stop condition */
           occnt_add,              /* flag for element name modification */
           b_indent[100],          /* list of indents */
           add_indent[];

extern double
           guess_orbit[6],
           al_errors[ALIGN_MAX],
           fd_errors[FIELD_MAX];

extern char
           line[MADX_LINE_MAX],
           tag_type[MAX_TAG][16],
           tag_code[MAX_TAG][16];

extern time_t last_time,
              start_time;

extern char filenames[100][500];
extern int  currentline[100];

extern double** trackstrarpositions; /* two dimensional array with track positions*/

#endif // MAD_STR_H


