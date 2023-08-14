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

ExpTree *derivative(ExpTree *expr, char *var) {
  if (expr == NULL) {
    return NULL;
  }

  switch (expr->type) {
    /* Derivative of a constant is 0 */
    case EXP_NUM:
      return newExpLeaf(EXP_NUM, "0");
    /* Derivative of a variable w.r.t itself is 1 */ 
    case EXP_VAR:
      if (strcmp(expr->data, var) == 0) {
        return newExpLeaf(EXP_NUM, "1");
    /* Derivative of a variable w.r.t another variable is 0 */  
      } else {
        return newExpLeaf(EXP_NUM, "0"); 
      }
    case EXP_ADD_OP:
      return newExpOp(EXP_ADD_OP, derivative(expr->left, var), derivative(expr->right, var));
    case EXP_SUB_OP:
      return newExpOp(EXP_SUB_OP, derivative(expr->left, var), derivative(expr->right, var));
    case EXP_MUL_OP:
      return newExpOp(EXP_ADD_OP,
        newExpOp(EXP_MUL_OP, derivative(expr->left, var), expr->right),
        newExpOp(EXP_MUL_OP, expr->left, derivative(expr->right, var)));
    case EXP_EXP_OP:
      return newExpOp(EXP_MUL_OP,
        newExpOp(EXP_MUL_OP, newExpLeaf(EXP_NUM, expr->right->data), derivative(expr->left, var)),
        newExpOp(EXP_EXP_OP, expr->left, newExpOp(EXP_SUB_OP, expr->right, newExpLeaf(EXP_NUM, "1"))));
    /* more cases for other operators */
    default:
      assert(false);
      return NULL;
  }
}







