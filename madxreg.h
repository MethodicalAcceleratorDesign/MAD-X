struct r_char_array             /* dynamic array of char for regex */
{
  int  max,                     /* max. array size */
       curr;                    /* current occupation */
  char* chars;
};

struct reg_token
{
  int type;                    /* 0: empty, 1: char, 2: simple,
                                  3: list,  4: dot */
  int rep;                     /* 0 no, 1 yes (including 0 repetitions) */
  int rep_cnt;                 /*  no. of repitions at current match */
  int rep_max;                 /* max. no. of repitions at current match */
  int invert;                  /* 0 no, 1 yes (for type 3 only) */
  int match_end;                /* 0 no, 1 yes (only last token for '$' */
  char c;                      /* for type = 1 */
  struct r_char_array* simple;
  struct r_char_array* list;
  struct reg_token* next;
  struct reg_token* previous;
};
