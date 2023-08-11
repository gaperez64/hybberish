#ifndef SYSLINEQ_H
#define SYSLINEQ_H

/* Ordinary differential expression trees */
typedef enum LINEQExpType {
  LINEQ_NUM,
  LINEQ_ADD_OP,
  LINEQ_SUB_OP,
  LINEQ_NEG,
  LINEQ_VAR,
} LINEQExpType;

typedef struct LINEQExpTree {
  char *data;
  LINEQExpType type;
  struct LINEQExpTree *left;
  struct LINEQExpTree *right;
} LINEQExpTree;

LINEQExpTree *newLineqExpLeaf(LINEQExpType, char *);
LINEQExpTree *newLineqExpOp(LINEQExpType, LINEQExpTree *, LINEQExpTree *);
LINEQExpTree *newLineqExpTree(LINEQExpType, char *, LINEQExpTree *,
                              LINEQExpTree *);
void delLineqExpTree(LINEQExpTree *);
void printLineqExpTree(LINEQExpTree *, FILE *);

/* Lists of LINEQs */
typedef struct LINEQList {
  LINEQExpTree *exp;
  struct LINEQList *next;
} LINEQList;

LINEQList *newLineqList(LINEQExpTree *);
/* Append the second one to the first list,
 * It assumes the second one is a single node!
 */
LINEQList *appLineqElem(LINEQList *, LINEQList *);
void delLineqList(LINEQList *);
void printLineqList(LINEQList *, FILE *);

#endif
