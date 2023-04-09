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
    char* data;
    ODEExpType type;
    struct ODEExpTree* left;
    struct ODEExpTree* right;
} ODEExpTree;
