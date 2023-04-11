#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "sysode.h"

ODEExpTree *newOdeExpLeaf(ODEExpType type, char *name) {
  ODEExpTree *tree = (ODEExpTree *)malloc(sizeof(ODEExpTree));
  tree->data = name;
  tree->type = type;
  switch(type) {
    case ODE_NUM:
      tree->left = NULL;
      tree->right = NULL;
      break;
    case ODE_VAR:
      tree->left = NULL;
      tree->right = NULL;
      break;
    default:
      assert(false);
      return NULL;
  }
  return tree;
}

ODEExpTree *newOdeExpOp(ODEExpType type,
                        ODEExpTree *left, ODEExpTree *right) {
  return newOdeExpTree(type, NULL, left, right);
}

ODEExpTree *newOdeExpTree(ODEExpType type, char *name, 
                          ODEExpTree *left, ODEExpTree *right) {
  ODEExpTree *tree = (ODEExpTree *)malloc(sizeof(ODEExpTree));
  tree->data = name;
  tree->type = type;
  switch(type) {
    /* binary operators */
    case ODE_ADD_OP:
    case ODE_SUB_OP:
    case ODE_MUL_OP:
    case ODE_DIV_OP:
    case ODE_EXP_OP:
      assert(left != NULL);
      tree->left = left;
      assert(right != NULL);
      tree->right = right;
      break;
    /* unary operators */
    case ODE_NEG:
    case ODE_FUN:
      assert(left != NULL);
      tree->left = left;
      assert(right == NULL);
      break;
    default:
      assert(false);
      return NULL;
  }
  return tree;
}

void delOdeExpTree(ODEExpTree *tree) {
  /* A simple depth-first search while freeing nodes post-order */
  if (tree->left != NULL)
    delOdeExpTree(tree->left);
  if (tree->right != NULL)
    delOdeExpTree(tree->right);
  /* Base case: delete the current node */
  if (tree->data != NULL)
    free(tree->data);
  free(tree);
}

static void printBinOp(ODEExpType type, FILE* where) {
  switch(type) {
    case ODE_ADD_OP:
      fprintf(where, " + ");
      break;
    case ODE_SUB_OP:
      fprintf(where, " - ");
      break;
    case ODE_MUL_OP:
      fprintf(where, " * ");
      break;
    case ODE_DIV_OP:
      fprintf(where, " / ");
      break;
    case ODE_EXP_OP:
      fprintf(where, " ^ ");
      break;
    default:
      assert(false);
  }
}

void printOdeExpTree(ODEExpTree *tree, FILE *where) {
  /* A simple depth-first search while printing in-order */
  switch(tree->type) {
    /* binary operators */
    case ODE_ADD_OP:
    case ODE_SUB_OP:
    case ODE_MUL_OP:
    case ODE_DIV_OP:
    case ODE_EXP_OP:
      fprintf(where, "(");
      assert(tree->left != NULL);
      printOdeExpTree(tree->left, where);
      printBinOp(tree->type, where);
      assert(tree->right != NULL);
      printOdeExpTree(tree->right, where);
      fprintf(where, ")");
      break;
    /* unary operators */
    case ODE_NEG:
      fprintf(where, "-");
      assert(tree->left != NULL);
      printOdeExpTree(tree->left, where);
      assert(tree->right == NULL);
      break;
    case ODE_FUN:
      assert(tree->data != NULL);
      fprintf(where, "%s(", tree->data);
      assert(tree->left != NULL);
      printOdeExpTree(tree->left, where);
      fprintf(where, ")");
      assert(tree->right == NULL);
      break;
    /* base cases */
    case ODE_NUM:
    case ODE_VAR:
      assert(tree->data != NULL);
      fprintf(where, "%s", tree->data);
      assert(tree->left == NULL);
      assert(tree->right == NULL);
      break;
    default:
      assert(false);
  }
}
