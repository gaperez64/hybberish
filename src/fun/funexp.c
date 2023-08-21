#include "funexp.h"
#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

ExpTree *newExpLeaf(ExpType type, char *name) {
  ExpTree *tree = (ExpTree *)malloc(sizeof(ExpTree));
  tree->data = name;
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
    copy->data = strdup(src->data);
	/* recursively copy the left tree node */
    copy->left = cpyExpTree(src->left);
	/* recursively copy the right tree node */
    copy->right = cpyExpTree(src->right);

    return copy;
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
        case EXP_ADD_OP:
            return newExpOp(EXP_ADD_OP, derivative(expr->left, var), derivative(expr->right, var));
        case EXP_SUB_OP:
            return newExpOp(EXP_SUB_OP, derivative(expr->left, var), derivative(expr->right, var));
        case EXP_MUL_OP:
            return newExpOp(EXP_ADD_OP,
                            newExpOp(EXP_MUL_OP, derivative(expr->left, var), cpyExpTree(expr->right)),
                            newExpOp(EXP_MUL_OP, cpyExpTree(expr->left), derivative(expr->right, var)));
        case EXP_EXP_OP:
            return newExpOp(EXP_MUL_OP,
                            newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, cpyExpTree(expr->right)->data), derivative(expr->left, var)),
                            newExpOp(EXP_EXP_OP, cpyExpTree(expr->left),
                                     newExpOp(EXP_SUB_OP, cpyExpTree(expr->right), newExpLeaf(EXP_NUM, "1"))));
        /* more cases for other operators */
        default:
            assert(false);
            return NULL;
    }
}

ExpTree *integral(ExpTree *expr, char *var) {
    if (expr == NULL) {
        return NULL;
    }
    
    switch (expr->type) {
        case EXP_NUM:
            return newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, expr->data), newExpLeaf(EXP_VAR, var));
        case EXP_VAR:
            if (strcmp(expr->data, var) == 0) {
                return newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, "0.5"), newExpOp(EXP_EXP_OP, newExpLeaf(EXP_VAR, var), newExpLeaf(EXP_NUM, "2")));
            } else {
                return newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, expr->data), newExpLeaf(EXP_VAR, var));
            }
        case EXP_ADD_OP:
            return newExpOp(EXP_ADD_OP, integral(expr->left, var), integral(expr->right, var));
        case EXP_SUB_OP: {
            ExpTree *left_integral = integral(expr->left, var);
            ExpTree *right_integral = integral(expr->right, var);
            return newExpOp(EXP_SUB_OP, left_integral, right_integral);
        }
        case EXP_MUL_OP: {
            ExpTree *left_integral = integral(expr->left, var);
            ExpTree *right = cpyExpTree(expr->right);
            return newExpOp(EXP_MUL_OP, left_integral, right);
        }
        case EXP_EXP_OP: {
            if (strcmp(expr->right->data, "1") == 0) {
                return integral(expr->left, var);
            } else {
                ExpTree *left = cpyExpTree(expr->left);
                ExpTree *right = cpyExpTree(expr->right);
                ExpTree *new_right = newExpOp(EXP_ADD_OP, right, newExpLeaf(EXP_NUM, "1"));
                ExpTree *power_rule = newExpOp(EXP_MUL_OP,
                                               newExpOp(EXP_DIV_OP, newExpLeaf(EXP_NUM, "1"), new_right),
                                               newExpOp(EXP_EXP_OP, left, new_right));
                 return power_rule;
            }
        }
        default:
            assert(false);
            return NULL;
    }
}






