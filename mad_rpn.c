#include "madx.h"

void
deco_init(void)
  /* initializes Polish decoding */
{
  expr_chunks = new_name_list("expr_chunks", 2000);
  cat = new_int_array(MAX_ITEM);
  deco = new_int_array(MAX_ITEM);
  d_var = new_int_array(MAX_ITEM);
  oper = new_int_array(MAX_ITEM);
  func = new_int_array(MAX_ITEM);
  cat_doubles = new_double_array(MAX_ITEM);
  doubles = new_double_array(MAX_D_ITEM);
  twiss_deltas = new_double_array(MAX_ITEM);
}

int
polish_expr(int c_item, char** item)   /* split input */
  /* deco output array containing:
     expression in Polish notation of length deco->curr,
     coded as 0-, 1+, 2*, 3/, 4^ (power),
     6 evaluate function
     100000000 + n = variable n (refers to vars),
     200000000 + n = function n (refers to functs),
     400000000 + n = real n (refers to doubles)
     -- Example: suppose a, b are variables 0 and 4, exp is function 3:
     then     3 * a * q[l] * q[k1] / exp((b - 1.57)^2) + 1.57
     would result in
     400000000 100000000 2 100000001 2 100000002 2
     100000003 400000001 0 400000002 3 200000003 3 400000001 1
     where 3 = real 0, 1.57 = real 1, 2 = real 2
     a = vars 0, q[l] vars 1, q[k1] vars 2, exp functs 3
  */
{
  int i, j, error, op, id, stack = 0, l_deco, l_double;
  int up[100][3] = {{-1, -1, -1}};

  l_deco = deco->curr = 0;
  l_double = doubles->curr;
  error = scan_expr(c_item, item);
  if (error) return error;
  for (i = 0; i < cat->curr; i++)
  {

    /* categories: 1: variable, 3: floating constant, 4: operator
       6: left par., 7: right par.     */
    switch (cat->i[i])
    {
      case 1:                              /* variable */
        if (l_deco == deco->max) grow_int_array(deco);
        deco->i[l_deco++] = 100000000 + d_var->i[i];
        break;
      case 3:                              /* constant */
        if (l_deco == deco->max) grow_int_array(deco);
        if (l_double == doubles->max) grow_double_array(doubles);
        doubles->a[l_double] = cat_doubles->a[i];
        deco->i[l_deco++] = 400000000 + l_double++;
        doubles->curr = l_double;
        break;
      case 4:
        if ((op = oper->i[i]) < 5)           /* operator */
        {
          id = op / 2;
          for (j = 2; j >= id; j--)
          {
            if (up[stack][j] > -1)
            {
              if (l_deco == deco->max) grow_int_array(deco);
              deco->i[l_deco++] = up[stack][j];
              up[stack][j] = -1;
            }
          }
          up[stack][id] = op;
        }
        else
        {
          if (l_deco == deco->max) grow_int_array(deco);
          deco->i[l_deco++] = 200000000 + func->i[i];  /* function */
        }
        break;
      case 6:      /*  '(' */
        stack++;
        for (j = 0; j < 3; j++)  up[stack][j] = -1;
        break;
      case 7:      /*  ')' */
        for (j = 2; j >= 0; j--)
        {
          if (up[stack][j] > -1)
          {
            if (l_deco == deco->max) grow_int_array(deco);
            deco->i[l_deco++] = up[stack][j];
          }
        }
        stack--;
        break;
      default:
        return 9;
    }   /* end switch */
  }     /* end loop over categories */
  for (j = 2; j >= 0; j--)   /* clear stack */
  {
    if (up[stack][j] > -1)
    {
      if (l_deco == deco->max) grow_int_array(deco);
      deco->i[l_deco++] = up[stack][j];
    }
  }
  deco->curr = l_deco;
  return 0;
}

double
polish_value(struct int_array* deco, char* expr_string)
  /* coded input (see below) */
  /* description see polish_expression */
{
  int i, k, kc, c_stack = -1;
  double stack[MAX_ITEM];
  char tmp[20];

  if (++polish_cnt > MAX_LOOP)
    fatal_error("circular call in expression", expr_string);
  stack[0] = 0;
  for (i = 0; i < deco->curr; i++)   /* decoding loop */
  {
    k = deco->i[i];
    if ( k < 5)     /* operator */
    {
      if (c_stack < 0)
      {
        fatal_error("stack underflow in expression:", expr_string);
      }
      else if (c_stack == 0)
      {
        stack[1] = stack[0]; stack[0] = 0;
      }
      else  c_stack--;


      switch(k)
      {
        case 0:
          stack[c_stack] -= stack[c_stack+1];
          break;
        case 1:
          stack[c_stack] += stack[c_stack+1];
          break;
        case 2:
          stack[c_stack] *= stack[c_stack+1];
          break;
        case 3:
          if (stack[c_stack+1] == 0.0)
          {
            warning("division by zero, result set to zero, expr:",
                    expr_string);
            stack[c_stack] = 0.0;
            break;
          }
          stack[c_stack] /= stack[c_stack+1];
          break;
        case 4:
          stack[c_stack] = pow(stack[c_stack],stack[c_stack+1]);
          break;
        default:
          fatal_error("illegal operator, Polish, expr.:", expr_string);
      }
    }
    else
    {
      kc = k / 100000000;  k -= 100000000 * kc;
      switch(kc)
      {
        case 1:            /* variable */
          stack[++c_stack] = act_value(k, expr_chunks);
          break;
        case 4:            /* real constant */
          stack[++c_stack] = doubles->a[k];
          break;
        case 2:            /* function */
          switch(k-1)      /* the offset is due to dummyfunction */
          {
            case 0:
              stack[c_stack] = fabs(stack[c_stack]);
              break;
            case 1:
              stack[c_stack] = sqrt(stack[c_stack]);
              break;
            case 2:
              stack[c_stack] = exp(stack[c_stack]);
              break;
            case 3:
              stack[c_stack] = log(stack[c_stack]);
              break;
            case 4:
              stack[c_stack] = log10(stack[c_stack]);
              break;
            case 5:
              stack[c_stack] = sin(stack[c_stack]);
              break;
            case 6:
              stack[c_stack] = cos(stack[c_stack]);
              break;
            case 7:
              stack[c_stack] = tan(stack[c_stack]);
              break;
            case 8:
              stack[c_stack] = asin(stack[c_stack]);
              break;
            case 9:
              stack[c_stack] = acos(stack[c_stack]);
              break;
            case 10:
              stack[c_stack] = atan(stack[c_stack]);
              break;
            case 11:
              stack[c_stack] = sinh(stack[c_stack]);
              break;
            case 12:
              stack[c_stack] = cosh(stack[c_stack]);
              break;
            case 13:
              stack[c_stack] = tanh(stack[c_stack]);
              break;
            case 14:
              stack[c_stack] = frndm();
              break;
            case 15:
              stack[c_stack] = grndm();
              break;
            case 16:
              stack[c_stack] = tgrndm(stack[c_stack]);
              break;
            case 17:
              stack[c_stack] = table_value();
              break;
            case 18: /* function "exist" */
              continue; /* value in stack not changed */
              break;
            case 19:
              stack[c_stack] = floor(stack[c_stack]);
              break;
            case 20:
              stack[c_stack] = ceil(stack[c_stack]);
              break;
            case 21:
              stack[c_stack] = rint(stack[c_stack]);
              break;
            case 22: {
              double int_part;
              stack[c_stack] = modf(stack[c_stack], &int_part);
              } break;
            default:
              fatal_error("polish_value: illegal function in expr:",
                          expr_string);
          }
          break;
        default:
          /* if you get here, then most likely someone has created
             more than 100000000 double precision numbers */
          sprintf(tmp, "%d", k-1);
          fatal_error("illegal type in Polish decoding: ", tmp);
          exit(1);
      }
    }
  }       /* end of decoding loop */
  polish_cnt--;
  return stack[0];
}

