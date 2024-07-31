#include "funexp.h"
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

ExpTree *newExpLeaf(ExpType type, char *name) {
  ExpTree *tree = (ExpTree *)malloc(sizeof(ExpTree));
  tree->data = name;
  // test code
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

ExpTree *newExpOp(ExpType type, ExpTree *left, ExpTree *right) {
  return newExpTree(type, NULL, left, right);
}

ExpTree *newExpTree(ExpType type, char *name, ExpTree *left, ExpTree *right) {
  ExpTree *tree = (ExpTree *)malloc(sizeof(ExpTree));
  tree->data = name;
  // tree->data = NULL;
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

void printExpTree(ExpTree *tree, FILE *where) {
  /* A simple depth-first search while printing in-order */
  assert(tree != NULL);
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

ExpTree *cpyExpTree(ExpTree *src) {
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

bool isLinear(ExpTree *expr) {
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

ExpTree *derivative(ExpTree *expr, char *var) {
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
    ExpTree *base = expr->left;
    ExpTree *exponent = expr->right;

    /* Derivative of base */
    ExpTree *base_derivative = derivative(base, var);

    /* Convert exponent from string to number */
    double exponent_val = atof(exponent->data);

    /* Derivative of exponent */
    ExpTree *exponent_minus_one =
        newExpOp(EXP_SUB_OP, exponent, newExpLeaf(EXP_NUM, "1"));
    ExpTree *exponent_derivative =
        newExpOp(EXP_MUL_OP, exponent, base_derivative);
    ExpTree *exponent_term = newExpOp(EXP_EXP_OP, base, exponent_minus_one);

    /* Convert (number - 1) back to a string */
    char exponent_str[50];
    snprintf(exponent_str, sizeof(exponent_str), "%.15g", exponent_val - 1);

    /* Replace original exponent with new string exponent */
    exponent_term->right = newExpLeaf(EXP_NUM, exponent_str);

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

    if (strcmp(function_name, "sin") == 0) {
      ExpTree *arg = expr->left;
      ExpTree *arg_derivative = derivative(arg, var);

      /* Derivative of sine */
      ExpTree *cos_func =
          newExpTree(EXP_FUN, strdup("cos"), cpyExpTree(arg), NULL);
      ExpTree *derivative_result = newExpOp(EXP_MUL_OP, cpyExpTree(cos_func),
                                            cpyExpTree(arg_derivative));

      return derivative_result;
    } else if (strcmp(function_name, "cos") == 0) {
      ExpTree *arg = expr->left;
      ExpTree *arg_derivative = derivative(arg, var);

      /* Derivative of cosine */
      ExpTree *neg_sin_func =
          newExpTree(EXP_FUN, strdup("sin"), cpyExpTree(arg), NULL);
      ExpTree *derivative_result =
          newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "-1"),
                   newExpOp(EXP_MUL_OP, cpyExpTree(neg_sin_func),
                            cpyExpTree(arg_derivative)));
      return derivative_result;
    } else if (strcmp(function_name, "sqrt") == 0) {
      ExpTree *arg = expr->left;
      ExpTree *arg_derivative = derivative(arg, var);

      /* Derivative of sqrt */
      ExpTree *half = newExpLeaf(EXP_NUM, "0.5");
      ExpTree *sqrt_derivative = newExpOp(
          EXP_MUL_OP, half,
          newExpOp(EXP_DIV_OP, cpyExpTree(arg_derivative),
                   newExpTree(EXP_FUN, strdup("sqrt"), cpyExpTree(arg), NULL)));

      return sqrt_derivative;
    }

  /* more cases for other operators */
  default:
    assert(false);
    return NULL;
  }
  }
}

ExpTree *integral(ExpTree *expr, char *var) {
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
    ExpTree *base = expr->left;
    ExpTree *exponent = expr->right;

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

    if (strcmp(function_name, "sin(x)") == 0) {
      /* Handle integral of sin(x) */
      ExpTree *integral_result = newExpOp(
          EXP_MUL_OP, newExpLeaf(EXP_NUM, "-1"),
          newExpTree(EXP_FUN, strdup("cos"), newExpLeaf(EXP_VAR, var), NULL));
      return integral_result;
    } else if (strstr(expr->data, "sin(") != NULL &&
               strstr(expr->data, "x)") != NULL) {
      char *n_str = strdup(expr->data + 4);
      n_str[strlen(n_str) - 1] = '\0';

      double n = atof(n_str);

      if (n != 0) {
        /* Handle integral of sin(nx) */
        char result_str[50];
        snprintf(result_str, sizeof(result_str), "(-1/%.1f)*cos(%.1f)", n, n);
        ExpTree *integral_result = newExpTree(EXP_FUN, strdup(result_str),
                                              newExpLeaf(EXP_VAR, var), NULL);
        free(n_str);
        return integral_result;
      }
    } else if (strcmp(function_name, "cos(x)") == 0) {
      /* Handle integral of cos(x) */
      ExpTree *integral_result = newExpOp(
          EXP_MUL_OP, newExpLeaf(EXP_NUM, "1"),
          newExpTree(EXP_FUN, strdup("sin"), newExpLeaf(EXP_VAR, var), NULL));
      return integral_result;
    } else if (strstr(expr->data, "cos(") != NULL &&
               strstr(expr->data, "x)") != NULL) {
      char *n_str = strdup(expr->data + 4);
      n_str[strlen(n_str) - 1] = '\0';

      double n = atof(n_str);

      if (n != 0) {
        /* Handle integral of cos(nx) */
        char result_str[50];
        snprintf(result_str, sizeof(result_str), "(1/%.1f)*sin(%.1f)", n, n);
        ExpTree *integral_result = newExpTree(EXP_FUN, strdup(result_str),
                                              newExpLeaf(EXP_VAR, var), NULL);
        free(n_str);
        return integral_result;
      }
    }

    else if (strcmp(function_name, "sqrt") == 0) {
      ExpTree *arg = expr->left;

      /* Compute the integral of sqrt(x) as (2/3) * x^(3/2) */
      ExpTree *integral_result = newExpOp(
          EXP_MUL_OP, newExpLeaf(EXP_NUM, "(2/3)"),
          newExpOp(EXP_EXP_OP, cpyExpTree(arg), newExpLeaf(EXP_NUM, "(3/2)")));

      return integral_result;
    }

  default:
    assert(false);
    return NULL;
  }
  }
}


bool isEqual(ExpTree *tree1, ExpTree *tree2) {
  /* Base case: empty trees are equal. */
  if (tree1 == NULL && tree2 == NULL)
    return true;

  /* Guard statements: trees not equal if current nodes not exactly equal. */
  // Both trees are (non-)empty
  if ((tree1 == NULL) != (tree2 == NULL))
    return false;
  if (tree1->type != tree2->type)
    return false;
  // Fail fast tactic: (not-)NULL-ness of left & right subtrees must match,
  //    though this is eventually caught by the base case as well
  if (((tree1->left == NULL) != (tree2->left == NULL)) ||
      ((tree1->right == NULL) != (tree2->right == NULL)))
    return false;
  if ((tree1->data == NULL) != (tree2->data == NULL))
    return false;
  else if (tree1->data != NULL && strcmp(tree1->data, tree2->data) != 0)
    return false;

  /* Recursive case: left & right subtrees respectively must be equal */
  return isEqual(tree1->left, tree2->left) &&
         isEqual(tree1->right, tree2->right);
}
