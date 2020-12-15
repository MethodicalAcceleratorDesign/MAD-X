#include "madx.h"

static int
v_length(const char* form)
{
  int ret = 0;
  if      (form[1] == 'I') sscanf(int_format, "%d", &ret);
  else if (form[1] == 'F') sscanf(float_format, "%d", &ret);
  return ret;
}

static int
get_val_num(char* in_string, int start, int end)
{
  int j, dot = 0, exp = 0, sign = 0;
  char c;
  for (j = start; j < end; j++)
  {
    c = in_string[j];
    if(!isdigit(c))
    {
      if ((c = in_string[j]) == '.')
      {
        if (dot || exp) return (j - 1);
        dot = 1;
      }
      else if (c == 'e' || c == 'd')
      {
        if (c == 'd') in_string[j] = 'e';
        if (exp) return (j - 1);
        else exp = j+1;
      }
      else if(strchr("+-", c))
      {
        if (exp != j || sign) return (j - 1);
        sign = 1;
      }
      else return (j - 1);
    }
  }
  return (j - 1);
}

// public interface

char*
v_format(const char* string)
  /* copies string to global variable var_form
     replacing  %S, %I, and %F by the user defined formats;
     %NF and %NI are replaced by the field lengths (!) of the defined formats */
{
  static char var_form[1000];

  const char *p, *q = string, *s = string, *t;
  char c;
  *var_form = '\0';
  while ((p = strpbrk(s, "NIFS")))
  {
    if (p > q)
    {
      t = p; t--;
      if (*t == '%')
      {
        c = *p;
        strncat(var_form, q, p - q);
        if (c == 'N')
        {
          sprintf(&var_form[strlen(var_form)], "%d", v_length(p));
          p++;
        }

        else if (c == 'F')  strcat(var_form, float_format);
        else if (c == 'S')  strcat(var_form, string_format);
        else if (c == 'I')  strcat(var_form, int_format);
        q = p; q++;
      }
    }
    s = ++p;
  }
  strcat(var_form, q);
  return var_form;
}

double
simple_double(char** toks, int start, int end)
{
  if (start > end && start + 1 != end)  return INVALID;
  else if (start == end) return myatof(toks[start]);
  else
  {
    if (*toks[start] == '-') return -myatof(toks[end]);
    else if (*toks[start] == '+') return myatof(toks[end]);
    else return INVALID;
  }
}
int
in_spec_list(char* string)
  /* checks for presence of special commands IF() etc. */
{
  char* cp;
  int i = 0, n = imin((int)strlen(string), 100);
  strncpy(c_dum->c, string, n); stolower(c_dum->c);
  supp_char(' ', c_dum->c);
  char* semicolon = strchr(c_dum->c, ';');
  while (special_comm_cnt[i])
  {
    if (special_comm_desc[i][0] == '>')
    {
      if ((cp = strchr(c_dum->c, special_comm_desc[i][1])) != NULL)
      {
        if (strncmp(++cp, &special_comm_desc[i][2], special_comm_cnt[i]) == 0
                && (semicolon == NULL || cp < semicolon))
            return i+1;
      }
    }
    else if (strncmp(c_dum->c, &special_comm_desc[i][0],special_comm_cnt[i])
             == 0)  return i+1;
    i++;
  }
  return 0;
}

char*
spec_join(char** it_list, int n)
  /* replaces variable in table(...) by original string */
{
  *c_join->c = '\0';
  if (n > 0) {
    assert(n < 1000);
    char* p[n];
    for (int j = 0; j < n; j++) p[j] = it_list[j];
    for (int j = 0; j < n; j++) {
      struct variable* var;
      if (strcmp(p[j], "table") == 0 && j+3 < n  && (var = find_variable(p[j+2], variable_list)) != NULL)
        p[j+2] = var->string;
    }
    for (int j = 0; j < n; j++) strcat(c_join->c, p[j]);
  }
  return c_join->c;
}

void
pre_split(char* inbuf, struct char_array* outbuf, int fill_flag)
  /* inserts blanks between tokens */
  /* fill_flag != 0 makes a 0 to be inserted into an empty "()" */
{
  char c, cp = ' ', cpnb = ' ', quote = ' ';
  int j, k, kn, sl = strlen(inbuf), cout = 0, quote_lv = 0, rb_level = 0;
  int left_b = 0, new_string = 1, c_digit = 0, f_equal = 0, comm_cnt = 0;
  int len = strlen(inbuf);

  while (2*len > outbuf->max) grow_char_array(outbuf);
  for (k = 0; k < sl; k++)
  {
    c = inbuf[k];
    if (quote_lv > 0)
    {
      if (c == quote)
      {
        quote_lv--; outbuf->c[cout++] = c; outbuf->c[cout++] = ' ';
      }
      else outbuf->c[cout++] = c == ' ' ? '@' : c;
    }
    else
    {
      c = inbuf[k];
      switch (c)
      {
        case '\"':
        case '\'':
          quote = c;
          quote_lv++; outbuf->c[cout++] = ' '; outbuf->c[cout++] = c;
          break;
        case '-':
          if (inbuf[k+1] == '>')
          {
            outbuf->c[cout++] = c; break;
          }
          /* FALLTHRU */
        case '+':
          if (left_b > 0)
          {
            outbuf->c[cout++] = ' ';
            outbuf->c[cout++] = '0';
            outbuf->c[cout++] = ' ';
            left_b = 0;
          }
          if (!(new_string > 0 && c_digit > 0 && strchr("ed",cp)) && cout > 0)
            outbuf->c[cout++] = ' ';
          outbuf->c[cout++] = c;
          if (!(new_string > 0 && c_digit > 0 && strchr("ed",cp)))
          {
            outbuf->c[cout++] = ' ';
            new_string = 1;
          }
          break;
        case '(':
          rb_level++;
          left_b = 1;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          outbuf->c[cout++] = c;
          outbuf->c[cout++] = ' ';
          break;
        case '>':
          if (cout > 0 && outbuf->c[cout-1] == '-')
          {
            outbuf->c[cout++] = c;
          }
          else
          {
            left_b = 0;
            new_string = 1;
            outbuf->c[cout++] = ' ';
            outbuf->c[cout++] = c;
            outbuf->c[cout++] = ' ';
          }
          break;
        case ')':
          rb_level--;
          if (fill_flag && cpnb == '(') outbuf->c[cout++] = '0';
          /* FALLTHRU */

        case '<':
        case ':':
        case '*':
        case '/':
        case '^':
        case '{':
        case '}':
        case '[':
        case ']':
        case '|':
        case '&':
          left_b = 0;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          outbuf->c[cout++] = c;
          outbuf->c[cout++] = ' ';
          break;
        case '=':
          f_equal = 1;
          left_b = 0;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          outbuf->c[cout++] = c;
          outbuf->c[cout++] = ' ';
          break;
        case ',': /* kept behind first "=", or if not first "," */
          /* not kept inside round brackets before '=' */
          left_b = 0;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          if (f_equal || (comm_cnt && rb_level == 0))
          {
            outbuf->c[cout++] = c;
            outbuf->c[cout++] = ' ';
          }
          comm_cnt++;
          break;
        case ';':
          left_b = 0;
          new_string = 1;
          outbuf->c[cout++] = ' ';
          break; 
        default:
          if (c == ' ') outbuf->c[cout++] = c;
          //This is to handle the case of hexdecimal numbers
          else if(c=='0' && inbuf[k+1]=='x' && k+1<sl && isalpha(inbuf[k-1])==0 && isdigit(inbuf[k-1])==0) {
            int start = k;
            outbuf->c[cout++] = ' ';
            for (int o=start; o < sl ; o++){
              if(inbuf[o]=='-' || inbuf[o]=='+'){
                outbuf->c[cout++] = inbuf[o];
                outbuf->c[cout++] = inbuf[o+1];
                if(isdigit(inbuf[o+2]) && o+2<sl) {
                  outbuf->c[cout++] = inbuf[o+2];
                  k = o+2;
                }
                else{
                  k = o+1;
                }
                outbuf->c[cout++] = ' ';
                
                break;
              }
              outbuf->c[cout++] = inbuf[o];
            }
          }
          else
          {
            left_b = 0;
            if (new_string && (isdigit(c) || c == '.'))
            {
              kn = get_val_num(inbuf, k, sl);
              for (j = k; j <= kn; j++) outbuf->c[cout++] = inbuf[j];
              outbuf->c[cout++] = ' ';
              k = kn;
            }
            else
            {
              new_string = 0;
              if (cout > 0 || c != ' ') outbuf->c[cout++] = c;
            }
          }
      }
      cp = c; if (c != ' ') cpnb = c;
    }
  }
  outbuf->c[cout] = '\0';
}


