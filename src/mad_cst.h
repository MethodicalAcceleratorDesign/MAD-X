#ifndef MAD_CST_H
#define MAD_CST_H

// types

struct expression;

struct constant
{
  char name[NAME_L];
  struct expression* exp;       /* pointer to defining expression (always) */
  double value;
  int stamp;
};

// interface

void get_defined_constants(void);

#endif // MAD_CST_H

