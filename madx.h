struct char_array           /* dynamic array of char */
{ 
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. array size */
       curr;                    /* current occupation */
  char* c;
};

struct char_array_list
{
  char name[NAME_L];
  int stamp;
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct char_array** ca;
};

struct char_p_array           /* dynamic array of char pointers */
{ 
  char name[NAME_L];
  int  max,                     /* max. array size */
       curr,                    /* current occupation */
       flag;                    /* ancillary flag */
  int stamp;
  char** p;
};

struct command                     /* holds one command */
{
  char name[NAME_L];
  char module[NAME_L];                 /* name of module it belongs to */
  char group[NAME_L];                  /* command group it belongs to */
  int stamp;
  int link_type;                       /* 0 none, 1 start, 2 end of group */
  int mad8_type;                       /* 0 none, else mad-8 element code */
  int beam_def;                        /* beam commands: 1 if defined */
  regex_t* reg_pattern;                /* for regular expressions */
  struct name_list* par_names;         /* names + input flag of parameters */
  struct command_parameter_list* par;  /* parameter pointer list */
};

struct command_list /* contains list of command pointers sorted by name */
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  int stamp;
  struct command** commands;    /* command pointer list */
};

struct command_list_list /* contains list of command lists */
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  struct command_list** command_lists;    /* command_list pointer list */
  int stamp;
};

struct command_parameter        /* holds one command parameter */
{
  char name[NAME_L];
  int type;                           /* 0 logical 1 integer 2 double 
                                         3 string 4 constraint */
                                      /* 11 int array 12 double array
                                         13 string array */
  int c_type;                         /* for type 4:
                                         1 min, 2 max, 3 both, 4 value */
  double double_value;                /* type 0, 1, 2, 4 */
  double c_min;                       /* type 4 */
  double c_max;                       /* type 4 */
  struct expression* expr;            /* type 1, 2, 4 */
  struct expression* min_expr;        /* type 4 */
  struct expression* max_expr;        /* type 4 */
  char* string;                       /* type 3 */
  int stamp;
  struct double_array* double_array;  /* type 11, 12 */
  struct expr_list* expr_list;        /* type 11, 12 */
  struct char_p_array* m_string;      /* type 13 */
  struct command_parameter* call_def; /* contains definitions for "bare"
                                         parameter input, e.g. option,echo */
};

struct command_parameter_list /* contains list of command parameter pointers */
{
  int stamp;
  char name[NAME_L];
  int  max,                             /* max. pointer array size */
       curr;                            /* current occupation */
  struct command_parameter** parameters;  /* command_parameter pointer list */
};

struct constant
{
  char name[NAME_L];
  struct expression* exp;       /* pointer to defining expression (always) */
  double value;
  int stamp;
};

struct constraint /* contains one constraint */
{
  char name[NAME_L];
  int  type;                    /* 1 minimum */
                                /* 2 maximum */
                                /* 3 both 1 + 2 */
                                /* 4 value */
  double value,
         c_min,
         c_max,
         weight;
  int stamp;
};

struct constraint_list /* contains list of constraints */
{
  int stamp;
  char name[NAME_L];
  int  max,                           /* max. pointer array size */
       curr;                          /* current occupation */
  struct constraint** constraints;    /* command pointer list */
};

struct double_array        /* dynamic array of double */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. array size */
       curr;                    /* current occupation */
  double* a;
};

struct element             /* each element is unique */
{
  char name[NAME_L];
  int def_type;                 /* 0 if defined separately, 
                                   1 if inside sequence */
  int bv;                       /* bv: 0 false, 1 true (invert angle for
                                   sequence bv = -1) */
  double length;
  struct command* def;          /* pointer to defining command */
  struct element* parent;       /* pointer to parent of element */
                                /* *this for base_type elements (rbend etc.) */
  int stamp;
  struct element* base_type;    /* pointer to base_type of element */
                                /* *this for base_type elements (rbend etc.) */
};

struct el_list /* contains list of element pointers sorted by name */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  struct element** elem;        /* element pointer list */
};

struct expression
{
  char name[NAME_L];
  char* string;                 /* expression in string form */
  int status;                   /* status flag: 0 not evaluated
                                                1 evaluated */
  struct int_array* polish;     /* pointer to Polish notation, or NULL */
  double value;                 /* actual value */
  int stamp;
};

struct expr_list
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct expression** list;     /* expression pointer list */
};

struct in_buffer
{
  char name[NAME_L];
  int flag;                    /* flag for logical tests */
  struct char_array* c_a;
  int stamp;
};

struct in_buff_list
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation = call level */
  FILE** input_files;           /* input file pointers */
  int stamp;
  struct in_buffer** buffers;     /* in_buff pointer list */
};

struct in_cmd          /* contains information about classified command */
{
  char name[NAME_L];
  char* label;         /* pointer to label: if != NULL then buffer this */
  int type;            /*    0 command from list;
                             1 link start command from list;
                             2 link end command from list;
                             3 element definition;
                             4 variable definition; */
  int sub_type;        /* position in cmd_match_base */
  int stamp;
  int decl_start;      /* start of declarative part in tok_list */
  int clone_flag;      /* if zero, clone can be dropped after decoding */
  struct char_p_array* tok_list; /* contains pointers to tokens */
  struct command* cmd_def;       /* points to command definition */
  struct command* clone;         /* points to clone of command definition */
};

struct in_cmd_list /* contains list of in_cmd pointers sorted by label */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* labels;     /* index list of labels */
  struct in_cmd** in_cmds;      /* in_cmd pointer list */
};

struct int_array           /* dynamic array of int */
{ 
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. array size */
       curr;                    /* current occupation */
  int* i;
};

struct macro     /* stores one line or macro definition */
{
  char name[NAME_L];
  int n_formal;                 /* no. of formal parameters */
  struct char_p_array* formal;  /* list of formal parameters */
  struct char_p_array* tokens;  /* token pointers into body if split (line) */
  struct char_array* body;      /* contains all statements */
  int stamp;
};

struct macro_list
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;
  struct macro** macros;
};

struct name_list /* contains list of index sorted names plus int inform. */
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  int* index;                   /* index for alphabetic access */
  int* inform;                  /* array parallel to names with integer */
  int stamp;
  char** names;                 /* element names for sort */
};

struct node                /* the sequence is a linked list of nodes */
{
  char name[NAME_L];
  char* base_name;           /* basic type */
  struct node* previous;
  struct node* next;
  int share;               /* 0 normal, 1 if shared */
  int occ_cnt;             /* element occurrence count at node */
  int obs_point;           /* observation point number (tracking) */
  int sel_err;             /* error select flag */
  int con_cnt;             /* constraint counter */
  int enable;              /* flag for correctors and monitors: 0 off, 1 on */
  double position;         /* s position in sequence [m] */
  double at_value;
  double length;
  double dipole_bv;        /* +1 or -1 (if beam_bv AND element_bv) */
  double other_bv;         /* equal to beam_bv (+1 or -1) */
  double chkick;           /* calculated by orbit correction module */
  double cvkick;           /* calculated by orbit correction module */
  int stamp;
  struct expression* at_expr;
  char* from_name;
  struct element* p_elem;  /* pointer to element if any */
  struct sequence* p_sequ;  /* pointer to sequence if any */
  struct double_array* p_al_err; /* pointer to alignment error array */
  struct double_array* p_fd_err; /* pointer to field error array */
  struct command* savebeta; /* pointer to savebeta command if any */
  struct constraint_list* cl; /* pointer to constraint list during match */
  struct double_array* obs_orbit; /* for track observation point */
};

struct node_list /* contains list of node pointers sorted by name */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of node (!) names */
                                /* node_name = el_name:occ_cnt */
  struct node** nodes;          /* node pointer list */
};

struct sequence
{
  /* original sequence */
  char name[NAME_L];
  char* refpos;                 /* reference position for insertion */
  int ref_flag;                 /* -1 for exit, 0 for centre, 1 for entry */
  int share;                    /* 0 normal, 1 if shared */
  int nested;                   /* 0 flat, 1 if nested */
  int con_cnt;                  /* constraint counter */
  int stamp;
  double length;                /* length as in declaration */
  struct expression* l_expr;    /* length expression as in declaration */
  struct node* start;           /* first node in sequence */
  struct node* end;             /* last node in sequence */
  struct node_list* nodes;      /* alphabetic list of nodes */
  struct el_list* cavities;     /* alphabetic list of cavities */
  struct command* beam;         /* pointer to beam attached */
  /* expanded sequence */
  int n_nodes;                  /* number of nodes when expanded */
  struct node* ex_start;        /* first node in expanded sequence */
  struct node* ex_end;          /* last node in expanded sequence */
  struct node* range_start;     /* first node of current range in sequence */
  struct node* range_end;       /* last node of current range in sequence */
  struct node_list* ex_nodes;   /* alphabetic list of nodes */
  struct table* tw_table;       /* pointer to latest twiss table created */
  struct constraint_list* cl;   /* pointer to constraint list during match */
};

struct sequence_list /* contains list of sequence pointers sorted by name */
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  struct sequence** sequs;      /* sequence pointer list */
  int stamp;
};

struct table
{
  char name[NAME_L];
  int  max,                     /* max. # rows */
       curr,                    /* current # rows */
       num_cols,                /* # columns - fixed */
       origin;                  /* 0 if created in job, 1 if read */
  struct char_p_array* header;  /* extra lines for file header */
  struct int_array* col_out;    /* column no.s to be written (in this order) */
  struct int_array* row_out;    /* flag for row: 1 write, 0 don't */
  struct char_p_array* node_nm; /* names of nodes at each row */
  struct char_p_array** l_head; /* extra lines to be put in front of a line */
  struct node** p_nodes;        /* pointers to nodes at each row */
  char*** s_cols;               /* string columns */
  double** d_cols;              /* double precision columns */
  int stamp;
  struct name_list* columns;    /* names + types (in inform): 
                                   1 double, 3 string */
  struct sequence* org_sequ;    /* pointer to sequence it refers to */
};

struct table_list
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* names;      /* index list of tables */
  struct table** tables;
  int stamp;
};

struct variable
{
  char name[NAME_L];
  int status;                 /* 0 value not evaluated, 1 evaluated */
  int type;                   /* 0 constant, 1 direct, 2 deferred, 3 string */
  int val_type;               /* 0 int 1 double (0..2) */
  char* string;               /* pointer to string if 3 */
  struct expression* expr;    /* pointer to defining expression (0..2) */
  double value;               /* (0..2) */
  int stamp;
};

struct var_list /* contains list of variable pointers sorted by name */
{
  int stamp;
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct name_list* list;       /* index list of names */
  struct variable** vars;       /* variable pointer list */
};

  struct val_mic {
         double before[2];
         double after[2];
  };

  struct id_mic {
         int   id_ttb;
         int   enable;
         struct val_mic val;
         struct node* p_node;
         struct id_mic *next;
         struct id_mic *previous;
  };

  struct orb_cor {
         double qx0;
         double qy0;
         struct id_mic *cor_table;
         struct id_mic *mon_table;
  };
