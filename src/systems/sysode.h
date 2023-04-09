enum ODEExpType {
  NUM,
  ADD_OP,
  SUB_OP,
  MUL_OP,
  DIV_OP,
  EXP_OP,
  NEG,
  VAR,
  FUN,
};

typedef struct ODEExpTree {
  char *data;
  ODEExpType type;
  struct ODEExpTree *left;
  struct ODEExpTree *right;
} ODEExpTree;

typedef struct ODEList {
  char *var;
  ODEExpTree *exp;
  struct ODEList *next;
} ODEList;

ODEExpTree *newOdeExpNode(char *, ODEExpType);
ODEExpTree *newOdeExpTree(char *, ODEExpType, ODEExpTree *, ODEExpTree *);
void delOdeExpTree(ODEExpTree *);
