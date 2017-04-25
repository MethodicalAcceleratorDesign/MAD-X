#ifndef MAD_TABLE_H
#define MAD_TABLE_H

// types

struct int_array;
struct char_p_array;

struct node;
struct in_cmd;
struct sequence;
struct name_list;
struct command_list;

struct table
{
  char name[NAME_L],
       type[NAME_L];            /* like "twiss", "survey" etc. */
  int  max,                     /* max. # rows */
       curr,                    /* current # rows */
       num_cols,                /* total # columns - fixed */
       org_cols,                /* original # columns from definition */
       dynamic,                 /* if != 0, values taken from current row */
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

struct table_list_list
{
  char name[NAME_L];
  int  max,                     /* max. pointer array size */
       curr;                    /* current occupation */
  struct table_list** table_lists;
  int stamp;
};

// interface

struct table*           make_table(const char* name, const char* type, const char* const *table_cols, const int* table_types, int rows);
struct table*           new_table(const char* name, const char* type, int rows, struct name_list* cols);
struct table_list*      new_table_list(int size);
struct table_list_list* new_table_list_list(int size);
struct table*           delete_table(struct table*);
struct table*           read_table(struct in_cmd*);
struct table*           read_my_table(struct in_cmd*);

void    check_table(char* string);
void    check_tabindex(char* string);
void    check_tabstring(char* string);
double  table_value(void);
void    table_add_header(struct table*, const char* format, ...);
void    add_to_table_list(struct table*, struct table_list*);
void    add_vars_to_table(struct table*, double scale);
void    set_vars_from_table(struct table*);
void    double_table(char* table);
void    grow_table(struct table*); /* doubles number of rows */
void    print_table(struct table*);
void    make_map_table(int* map_table_max_rows);
int     get_table_range(char* range, struct table*, int* rows);
void    out_table(const char* tname, struct table*, const char* filename);
void    reset_count(const char* table); /* resets table counter to zero */
void    sector_out(char* sector_table_name, double* pos, double* kick, double* rmatrix, double* tmatrix);
void    table_range(char* table, char* range, int* rows);

void    rename_table(struct table *tbl, const char *name );
int     remove_table_from_table_list(const char *name, struct table_list* tl);
struct table *detach_table_from_table_list(const char *name, struct table_list* tl);

void    augment_count(const char* table);
void    augmentcountonly(const char* table);

//int     str_from_table     (const char* table, const char* name, int* row, char* val);
//int     str_from_tablet    (struct table *tbl, const char* name, int* row, char* val);
//int     nodename_from_table_row(const char* table, /* no name   */ const int* row, char* string);

int     table_length(const char* table);
int     table_exists(const char* table);
int     table_column_exists(const char* table, const char* name);
int     table_cell_exists(const char* table, const char* name, const int* row);
int     table_header_exists(const char* table, const char *name);

int     double_from_table_header(const char* table, const char* name, double* val);

int     double_from_table_row(const char* table, const char* name, const int* row, double* val);
int     string_from_table_row(const char* table, const char* name, const int* row, char* string);

int     double_to_table_row  (const char* table, const char* name, const int* row, const double* val);
int     string_to_table_row  (const char* table, const char* name, const int* row, const char* string);

int     double_to_table_curr (const char* table, const char* name, const double* val);
int     vector_to_table_curr (const char* table, const char* name, const double* vals, const int* nval);
int     string_to_table_curr (const char* table, const char* name, const char* string);
int     comment_to_table_curr(const char* table, const char* comment, const int* length);

// double  get_table_value(const char* table_s, const char *row_s, const char *col_s); // not used
// void    set_table_value(const char* table_s, const char *row_s, const char *col_s, double *val); // not used

struct column_info{
  void* data;
  int   length;
  char  datatype;
  char  datasize;
};

struct column_info   table_get_column(char* table_name,char* column_name);
struct char_p_array *table_get_header(char* table_name);

#endif // MAD_TABLE_H

