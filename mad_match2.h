#ifndef MAD_MATCH2_H
#define MAD_MATCH2_H

// types

struct node;
struct in_cmd;
struct expression;

// interface

void match2_match(struct in_cmd*);
void match2_end(struct in_cmd*);
void match2_macro(struct in_cmd*);
void match2_constraint(struct in_cmd*);
void match2_setconstrinrange(struct node**, double, char*, char, char*);
void match2_disasambleconstraint(struct in_cmd*);
int  match2_evaluate_exressions(int, int, double*);
void match2_delete_expressions(void);
int  match2_augmentnmacros(void);      /*increases space in the working arrays*/
int  match2_augmentnconstraints(void); /*increases space in the working arrays*/

/* frs 10.03.2008 */
int  match2_print_var(struct in_cmd*);
void match2_delete_arrays();
void match2_alloc_arrays();
void match2_init_arrays();

// variables

extern int MAX_MATCH_CONS; /*these are set to proper values at the initialization of the match2 module*/
extern int MAX_MATCH_MACRO; /*zero values means that it is not initialized yet*/

extern char match2_keepexpressions; /*do not delete expressions at the end matching used by match with PTC knobs*/

extern char   **match2_macro_name;
extern char*  **match2_cons_name;
extern double **match2_cons_value;
extern double **match2_cons_value_rhs;
extern double **match2_cons_value_lhs;
extern double **match2_cons_weight;
extern char   **match2_cons_sign;
extern struct expression* **match2_cons_rhs;
extern struct expression* **match2_cons_lhs;
extern int      match2_cons_curr[3];

#endif // MAD_MATCH2_H


