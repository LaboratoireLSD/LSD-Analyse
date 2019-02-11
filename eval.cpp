#include "eval.h"

char delims[] = "+-*/^)(";

int eval::evaluate(char *line, double *val)
{
      double arg;
      char *ptr = line, *str, *endptr;
      int ercode;
      
      state = 0;

      pack(line);

      while (*ptr)
      {
            switch (state)
            {
            case 0:
                  if (NULL != (str = getexp(ptr)))
                  {
                        if ('(' == *str)
                        {
                              push_op(*str);
                              ptr += strlen(str);
                              break;
                        }

                        if (0.0 == (arg = strtod(str, &endptr)) &&
                              NULL == strchr(str, '0'))
                        {
                              return ERROR;
                        }
                        push_arg(arg);
                        ptr += strlen(str);
                  }
                  else  return ERROR;

                  state = 1;
                  break;

            case 1:
                  if (NULL == (str = getop(ptr)))
                        return ERROR;

                  if (strchr(delims, *str))
                  {
                        if (')' == *str)
                        {
                              if (SUCCESS > (ercode = do_paren()))
                                    return ercode;
                        }
                        else
                        {
                              push_op(*str);
                              state = 0;
                        }

                        ptr += strlen(str);
                  }
                  else  return ERROR;

                  break;
            }
      }

      while (1 < arg_sptr)
      {
            if (SUCCESS > (ercode = do_op()))
                  return ercode;
      }
      if (!op_sptr)
            return pop_arg(val);
      else  return ERROR;
}

/*
**  Evaluate stacked arguments and operands
*/

int eval::do_op(void)
{
      double arg1, arg2;
      int op;

      if (ERROR == pop_op(&op))
            return ERROR;

      pop_arg(&arg1);
      pop_arg(&arg2);

      switch (op)
      {
      case '+':
            push_arg(arg2 + arg1);
            break;

      case '-':
            push_arg(arg2 - arg1);
            break;

      case '*':
            push_arg(arg2 * arg1);
            break;

      case '/':
            if (0.0 == arg1)
                  return R_ERROR;
            push_arg(arg2 / arg1);
            break;

      case '^':
            if (0.0 > arg2)
                  return R_ERROR;
            push_arg(pow(arg2, arg1));
            break;

      case '(':
            arg_sptr += 2;
            break;

      default:
            return ERROR;
      }
      if (1 > arg_sptr)
            return ERROR;
      else  return op;
}

/*
**  Evaluate one level
*/

int eval::do_paren(void)
{
      int op;

      if (1 > parens--)
            return ERROR;
      do
      {
            if (SUCCESS > (op = do_op()))
                  break;
      } while ('('!= op);
      return op;
}

/*
**  Stack operations
*/

void eval::push_op(char op)
{
      if ('(' == op)
            ++parens;
      op_stack[op_sptr++] = op;
}

void eval::push_arg(double arg)
{
      arg_stack[arg_sptr++] = arg;
}

STATUS eval::pop_arg(double *arg)
{
      *arg = arg_stack[--arg_sptr];
      if (0 > arg_sptr)
            return ERROR;
      else  return SUCCESS;
}

STATUS eval::pop_op(int *op)
{
      if (!op_sptr)
            return ERROR;
      *op = op_stack[--op_sptr];
      return SUCCESS;
}

/*
**  Get an expression
*/

char* eval::getexp(char *str)
{
      char *ptr = str, *tptr = token;

      while (*ptr)
      {
            if (strchr(delims, *ptr))
            {
                  if ('-' == *ptr)
                  {
                        if (str != ptr && 'E' != ptr[-1])
                              break;
                  }

                  else if (str == ptr)
                        return getop(str);

                  else if ('E' == *ptr)
                  {
                        if (!isdigit(ptr[1]) && '-' != ptr[1])
                              return NULL;
                  }
                  else break;
            }

            *tptr++ = *ptr++;
      }
      *tptr = NUL;

      return token;
}

/*
**  Get an operator
*/

char* eval::getop(char *str)
{
      *token = *str;
      token[1] = NUL;
      return token;
}

/*
**  Remove whitespace & capitalize
*/

void eval::pack(char *str)
{
      char *ptr = str, *p;

      strupr(str);

      for ( ; *ptr; ++ptr)
      {
            p = ptr;
            while (*p && isspace(*p))
                  ++p;
            if (ptr != p)
                  strcpy(ptr, p);
      }
}

void eval::strupr(char *s)
{
  while (*s)
    {
      *s = toupper((unsigned char) *s);
      s++;
    }
}
