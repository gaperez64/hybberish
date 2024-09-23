#include "funexp.h"
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

ExpTree *newExpLeaf(const ExpType type, const char *const name) {
  /* Leaves are fully dependent on their data, so it must be specified. */
  assert(name != NULL);

  ExpTree *tree = (ExpTree *)malloc(sizeof(ExpTree));
  tree->data = strdup(name);
  tree->type = type;
  switch (type) {
  case EXP_NUM:
    tree->left = NULL;
    tree->right = NULL;
    break;
  case EXP_VAR:
    tree->left = NULL;
    tree->right = NULL;
    break;
  default:
    assert(false);
    return NULL;
  }
  return tree;
}

ExpTree *newExpOp(const ExpType type, ExpTree *left, ExpTree *right) {
  return newExpTree(type, NULL, left, right);
}

ExpTree *newExpTree(const ExpType type, char *name, ExpTree *left,
                    ExpTree *right) {
  ExpTree *tree = (ExpTree *)malloc(sizeof(ExpTree));
  tree->data = name;
  tree->type = type;
  switch (type) {
  /* binary operators */
  case EXP_ADD_OP:
  case EXP_SUB_OP:
  case EXP_MUL_OP:
  case EXP_DIV_OP:
  case EXP_EXP_OP:
    assert(left != NULL);
    tree->left = left;
    assert(right != NULL);
    tree->right = right;
    break;
  /* unary operators */
  case EXP_NEG:
  case EXP_FUN:
    assert(left != NULL);
    tree->left = left;
    assert(right == NULL);
    tree->right = right;
    break;
  default:
    assert(false);
    return NULL;
  }
  return tree;
}

void delExpTree(ExpTree *tree) {
  assert(tree != NULL);

  /* A simple depth-first search while freeing nodes post-order */
  if (tree->left != NULL)
    delExpTree(tree->left);
  if (tree->right != NULL)
    delExpTree(tree->right);
  /* Base case: delete the current node */
  if (tree->data != NULL)
    free(tree->data);
  free(tree);
}

static void printBinOp(ExpType type, FILE *where) {
  switch (type) {
  case EXP_ADD_OP:
    fprintf(where, " + ");
    break;
  case EXP_SUB_OP:
    fprintf(where, " - ");
    break;
  case EXP_MUL_OP:
    fprintf(where, " * ");
    break;
  case EXP_DIV_OP:
    fprintf(where, " / ");
    break;
  case EXP_EXP_OP:
    fprintf(where, "^");
    break;
  default:
    assert(false);
  }
}

void printExpTree(const ExpTree *tree, FILE *where) {
  /* A simple depth-first search while printing in-order */
  assert(tree != NULL);
  assert(where != NULL);

  switch (tree->type) {
  /* binary operators */
  case EXP_ADD_OP:
  case EXP_SUB_OP:
  case EXP_MUL_OP:
  case EXP_DIV_OP:
  case EXP_EXP_OP:
    fprintf(where, "(");
    assert(tree->left != NULL);
    printExpTree(tree->left, where);
    printBinOp(tree->type, where);
    assert(tree->right != NULL);
    printExpTree(tree->right, where);
    fprintf(where, ")");
    break;
  /* unary operators */
  case EXP_NEG:
    fprintf(where, "-");
    assert(tree->left != NULL);
    printExpTree(tree->left, where);
    assert(tree->right == NULL);
    break;
  case EXP_FUN:
    assert(tree->data != NULL);
    fprintf(where, "%s(", tree->data);
    assert(tree->left != NULL);
    printExpTree(tree->left, where);
    fprintf(where, ")");
    assert(tree->right == NULL);
    break;
  /* base cases */
  case EXP_NUM:
  case EXP_VAR:
    assert(tree->data != NULL);
    fprintf(where, "%s", tree->data);
    assert(tree->left == NULL);
    assert(tree->right == NULL);
    break;
  default:
    assert(false);
  }
}

ExpTree *cpyExpTree(const ExpTree *const src) {
  if (src == NULL) {
    return NULL;
  }

  ExpTree *copy = (ExpTree *)malloc(sizeof(ExpTree));
  if (copy == NULL) {
    fprintf(stderr, "Memory allocation error\n");
    exit(EXIT_FAILURE);
  }

  copy->type = src->type;
  copy->data = (src->data != NULL) ? strdup(src->data) : NULL;

  /* recursively copy the left tree node */
  copy->left = cpyExpTree(src->left);
  /* recursively copy the right tree node */
  copy->right = cpyExpTree(src->right);

  return copy;
}

bool isLinear(const ExpTree *expr) {
  if (expr == NULL) {
    return false;
  }

  switch (expr->type) {
  case EXP_VAR:
    return true;

  case EXP_MUL_OP:
    /* Check if it's of the form 'n * x' where n is a constant */
    if (expr->left != NULL && expr->left->type == EXP_NUM &&
        expr->right != NULL && expr->right->type == EXP_VAR) {
      return true;
    }
    break;

  case EXP_ADD_OP:
    /* Check if both sides of the addition are linear expressions */
    return isLinear(expr->left) && isLinear(expr->right);

  case EXP_FUN:
    /* Check for sin and cos functions */
    if (strcmp(expr->data, "sin") == 0 || strcmp(expr->data, "cos") == 0) {
      return isLinear(expr->left);
    }
    break;

  default:
    break;
  }

  return false;
}

ExpTree *derivative(const ExpTree *expr, const char *var) {
  if (expr == NULL) {
    return NULL;
  }

  switch (expr->type) {
  case EXP_NUM:
    return newExpLeaf(EXP_NUM, "0");
  case EXP_VAR:
    if (strcmp(expr->data, var) == 0) {
      return newExpLeaf(EXP_NUM, "1");
    } else {
      return newExpLeaf(EXP_NUM, "0");
    }
  case EXP_ADD_OP: {
    ExpTree *left_derivative = derivative(expr->left, var);
    ExpTree *right_derivative = derivative(expr->right, var);
    return newExpOp(EXP_ADD_OP, left_derivative, right_derivative);
  }
  case EXP_SUB_OP: {
    ExpTree *left_derivative = derivative(expr->left, var);
    ExpTree *right_derivative = derivative(expr->right, var);
    return newExpOp(EXP_SUB_OP, left_derivative, right_derivative);
  }
  case EXP_MUL_OP: {
    ExpTree *left_derivative = derivative(expr->left, var);
    ExpTree *right_derivative = derivative(expr->right, var);
    ExpTree *left_copy = cpyExpTree(expr->left);
    ExpTree *right_copy = cpyExpTree(expr->right);

    ExpTree *left_term = newExpOp(EXP_MUL_OP, left_derivative, right_copy);
    ExpTree *right_term = newExpOp(EXP_MUL_OP, left_copy, right_derivative);

    return newExpOp(EXP_ADD_OP, left_term, right_term);
  }

  case EXP_EXP_OP: {
    ExpTree *base = cpyExpTree(expr->left);
    ExpTree *exponent = cpyExpTree(expr->right);

    /* Derivative of base */
    ExpTree *base_derivative = derivative(base, var);

    /* Convert exponent from string to number */
    double exponent_val = atof(exponent->data);

    /* Convert (number - 1) back to a string */
    char exponent_str[50];
    snprintf(exponent_str, sizeof(exponent_str), "%.15g", exponent_val - 1);

    /* Derivative of exponent */
    ExpTree *exponent_derivative =
        newExpOp(EXP_MUL_OP, exponent, base_derivative);
    ExpTree *exponent_term =
        newExpOp(EXP_EXP_OP, base, newExpLeaf(EXP_NUM, exponent_str));

    /* Final derivative result */
    ExpTree *final_derivative =
        newExpOp(EXP_MUL_OP, exponent_derivative, exponent_term);

    return final_derivative;
  }

  case EXP_FUN: {
    char *function_name = strdup(expr->data);
    for (int i = 0; function_name[i]; i++) {
      function_name[i] = tolower(function_name[i]);
    }

    ExpTree *derivative_result = NULL;

    if (strcmp(function_name, "sin") == 0) {
      ExpTree *arg = expr->left;
      ExpTree *arg_derivative = derivative(arg, var);

      /* Derivative of sine */
      ExpTree *cos_func =
          newExpTree(EXP_FUN, strdup("cos"), cpyExpTree(arg), NULL);
      derivative_result = newExpOp(EXP_MUL_OP, cos_func, arg_derivative);
    } else if (strcmp(function_name, "cos") == 0) {
      ExpTree *arg = expr->left;
      ExpTree *arg_derivative = derivative(arg, var);

      /* Derivative of cosine */
      ExpTree *neg_sin_func =
          newExpTree(EXP_FUN, strdup("sin"), cpyExpTree(arg), NULL);
      derivative_result =
          newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "-1"),
                   newExpOp(EXP_MUL_OP, neg_sin_func, arg_derivative));
    } else if (strcmp(function_name, "sqrt") == 0) {
      ExpTree *arg = expr->left;
      ExpTree *arg_derivative = derivative(arg, var);

      /* Derivative of sqrt */
      ExpTree *half = newExpLeaf(EXP_NUM, "0.5");
      derivative_result = newExpOp(
          EXP_MUL_OP, half,
          newExpOp(EXP_DIV_OP, arg_derivative,
                   newExpTree(EXP_FUN, strdup("sqrt"), cpyExpTree(arg), NULL)));
    }

    /* Clean. */
    free(function_name);

    return derivative_result;
  }

  /* more cases for other operators */
  default:
    assert(false);
    return NULL;
  }
}

ExpTree *integral(const ExpTree *expr, const char *var) {
  if (expr == NULL) {
    return NULL;
  }

  switch (expr->type) {
  case EXP_NUM:
    return newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, expr->data),
                    newExpLeaf(EXP_VAR, var));
  case EXP_VAR:
    if (strcmp(expr->data, var) == 0) {
      return newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "0.5"),
                      newExpOp(EXP_EXP_OP, newExpLeaf(EXP_VAR, var),
                               newExpLeaf(EXP_NUM, "2")));
    } else {
      return newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, expr->data),
                      newExpLeaf(EXP_VAR, var));
    }
  case EXP_ADD_OP:
    return newExpOp(EXP_ADD_OP, integral(expr->left, var),
                    integral(expr->right, var));

  case EXP_SUB_OP: {
    ExpTree *left_integral = integral(expr->left, var);
    ExpTree *right_integral = integral(expr->right, var);
    return newExpOp(EXP_SUB_OP, left_integral, right_integral);
  }

  case EXP_MUL_OP: {
    ExpTree *left_integral = integral(expr->left, var);
    ExpTree *right = integral(expr->right, var);
    return newExpOp(EXP_MUL_OP, left_integral, right);
  }

  case EXP_EXP_OP: {
    ExpTree *base = cpyExpTree(expr->left);
    ExpTree *exponent = cpyExpTree(expr->right);

    /* Convert exponent from string to number */
    double exponent_val = atof(exponent->data);

    /* Integral of exponent */
    ExpTree *exponent_plus_one =
        newExpOp(EXP_ADD_OP, exponent, newExpLeaf(EXP_NUM, "1"));
    ExpTree *exponent_integral =
        newExpOp(EXP_DIV_OP, newExpLeaf(EXP_NUM, "1"), exponent_plus_one);
    ExpTree *exponent_term = newExpOp(EXP_EXP_OP, base, exponent_plus_one);

    /* Convert (number + 1) back to a string */
    char exponent_str[50];
    snprintf(exponent_str, sizeof(exponent_str), "%.15g", exponent_val + 1);

    /* Replace original exponent with new string exponent */
    exponent_term->right = newExpLeaf(EXP_NUM, exponent_str);

    /* Final integral result */
    ExpTree *final_integral =
        newExpOp(EXP_MUL_OP, exponent_integral, exponent_term);

    return final_integral;
  }

  case EXP_FUN: {
    char *function_name = strdup(expr->data);
    for (int i = 0; function_name[i]; i++) {
      function_name[i] = tolower(function_name[i]);
    }

    ExpTree *integral_result = NULL;

    if (strcmp(function_name, "sin(x)") == 0) {
      /* Handle integral of sin(x) */
      integral_result = newExpOp(
          EXP_MUL_OP, newExpLeaf(EXP_NUM, "-1"),
          newExpTree(EXP_FUN, strdup("cos"), newExpLeaf(EXP_VAR, var), NULL));
    } else if (strstr(expr->data, "sin(") != NULL &&
               strstr(expr->data, "x)") != NULL) {
      char *n_str = strdup(expr->data + 4);
      n_str[strlen(n_str) - 1] = '\0';

      double n = atof(n_str);

      if (n != 0) {
        /* Handle integral of sin(nx) */
        char result_str[50];
        snprintf(result_str, sizeof(result_str), "(-1/%.1f)*cos(%.1f)", n, n);
        integral_result = newExpTree(EXP_FUN, strdup(result_str),
                                     newExpLeaf(EXP_VAR, var), NULL);
        free(n_str);
      }
    } else if (strcmp(function_name, "cos(x)") == 0) {
      /* Handle integral of cos(x) */
      integral_result = newExpOp(
          EXP_MUL_OP, newExpLeaf(EXP_NUM, "1"),
          newExpTree(EXP_FUN, strdup("sin"), newExpLeaf(EXP_VAR, var), NULL));
    } else if (strstr(expr->data, "cos(") != NULL &&
               strstr(expr->data, "x)") != NULL) {
      char *n_str = strdup(expr->data + 4);
      n_str[strlen(n_str) - 1] = '\0';

      double n = atof(n_str);

      if (n != 0) {
        /* Handle integral of cos(nx) */
        char result_str[50];
        snprintf(result_str, sizeof(result_str), "(1/%.1f)*sin(%.1f)", n, n);
        integral_result = newExpTree(EXP_FUN, strdup(result_str),
                                     newExpLeaf(EXP_VAR, var), NULL);
        free(n_str);
      }
    }

    else if (strcmp(function_name, "sqrt") == 0) {
      ExpTree *arg = expr->left;

      /* Compute the integral of sqrt(x) as (2/3) * x^(3/2) */
      integral_result = newExpOp(
          EXP_MUL_OP, newExpLeaf(EXP_NUM, "(2/3)"),
          newExpOp(EXP_EXP_OP, cpyExpTree(arg), newExpLeaf(EXP_NUM, "(3/2)")));
    }

    /* Clean. */
    free(function_name);

    return integral_result;
  }

  default:
    assert(false);
    return NULL;
  }
}

ExpTree *definiteIntegral(const ExpTree *expr, const char *var,
                          const ExpTree *lowerBound,
                          const ExpTree *upperBound) {
  assert(expr != NULL);
  assert(var != NULL);
  assert(lowerBound != NULL);
  assert(upperBound != NULL);
  assert(lowerBound->type == EXP_NUM || lowerBound->type == EXP_VAR);
  assert(upperBound->type == EXP_NUM || upperBound->type == EXP_VAR);

  ExpTree *integrated = integral(expr, var);
  ExpTree *lowerSubstituted = substitute(integrated, var, lowerBound);
  ExpTree *upperSubstituted = substitute(integrated, var, upperBound);
  ExpTree *definite = newExpOp(EXP_SUB_OP, upperSubstituted, lowerSubstituted);

  /* Clean */
  delExpTree(integrated);

  return definite;
}

bool isEqual(const ExpTree *expr1, const ExpTree *expr2) {
  /* Base case: empty trees are equal. */
  if (expr1 == NULL && expr2 == NULL)
    return true;

  /* Guard statements: two trees not equal if the current nodes
    are not exactly equal/identical. */
  // Both trees must be (non-)empty.
  if ((expr1 == NULL) != (expr2 == NULL))
    return false;
  // The node types must match.
  if (expr1->type != expr2->type)
    return false;
  // Fail fast tactic: (not-)NULL-ness of left & right subtrees must match,
  //    though this is eventually caught by the base case as well
  if (((expr1->left == NULL) != (expr2->left == NULL)) ||
      ((expr1->right == NULL) != (expr2->right == NULL)))
    return false;
  // Both trees must have (non-)empty data.
  if ((expr1->data == NULL) != (expr2->data == NULL))
    return false;
  // The data of both trees must match.
  else if (expr1->data != NULL && strcmp(expr1->data, expr2->data) != 0)
    return false;

  /* Recursive case: left & right subtrees respectively must be equal */
  return isEqual(expr1->left, expr2->left) &&
         isEqual(expr1->right, expr2->right);
}

unsigned int degreeMonomial(const ExpTree *expr) {
  assert(expr != NULL);

  switch (expr->type) {
  // Base case: A number constant is of degree zero.
  case EXP_NUM:
    return 0;

  // Base case: A variable is of degree 1.
  case EXP_VAR:
    return 1;

  case EXP_NEG:
    assert(expr->left != NULL);
    assert(expr->right == NULL);
    return degreeMonomial(expr->left);

  case EXP_MUL_OP:
    assert(expr->left != NULL);
    assert(expr->right != NULL);
    return degreeMonomial(expr->left) + degreeMonomial(expr->right);

  case EXP_EXP_OP:
    assert(expr->left != NULL);
    assert(expr->right != NULL);
    /* Restrict exponents to non-negative integers. */
    assert(expr->right->type == EXP_NUM);
    assert(expr->left->type == EXP_VAR);

    const char *pattern = "^[0-9]*(\\.0*)?$";
    regex_t regex;
    assert(!regcomp(&regex, pattern, REG_EXTENDED));
    assert(!regexec(&regex, expr->right->data, 0, NULL, 0));
    regfree(&regex);

    /* Do rounding to avoid floating point errors, only integer degrees are
     * wanted. */
    double exponent = round(atof(expr->right->data));
    assert(exponent >= 0);

    // TODO: This does not make sense. This code fails an assertion
    // for (2^4). Also, is it possible for (x^3)^4 to be passed in?
    // If not, then specify a pre-condition to prevent it, e.g. that
    // the expression must be a monomial or smth?
    // Also, get rid of the regex; the parser restricts exponents to
    // integers.
    return (unsigned int)exponent;

  /* Invalid subexpression for a monomial. */
  default:
    assert(false);
    break;
  }
}
