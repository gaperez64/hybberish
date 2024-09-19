/**
 * @file funexp.h
 * @brief An interface for the creation and manipulation of expression trees.
 * @version 0.1
 * @date 2024-09-18
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef FUNEXP_H
#define FUNEXP_H

#include <math.h>
#include <regex.h>
#include <stdbool.h>
#include <stdio.h>

/**
 * @brief An enumeration of expression tree node types.
 */
typedef enum ExpType {
  EXP_NUM,    ///< A **leaf node** type: a number constant.
  EXP_VAR,    ///< A **leaf node** type: a variable.
  EXP_ADD_OP, ///< An **internal node** type: a binary add (a + b) operator.
  EXP_SUB_OP, ///< An **internal node** type: a binary sub (a - b) operator.
  EXP_MUL_OP, ///< An **internal node** type: a binary mul (a * b) operator.
  EXP_DIV_OP, ///< An **internal node** type: a binary div (a / b) operator.
  EXP_EXP_OP, ///< An **internal node** type: a binary exp (a ^ b) operator.
  EXP_NEG,    ///< An **internal node** type: a unary neg (- a) operator.
  EXP_FUN,    ///< An **internal node** type: an arbirary funtion.
} ExpType;

/**
 * @brief A binary expression tree node.
 * @details The node type has a large impact on the following aspects:
 *    - The NULL-ness or the valid values of the data member.
 *    - The NULL-ness or non-NULL-ness of the left and or right subtrees.
 *
 * Well-formedness of an expression tree node is only guaranteed if the
 * correct constructors are used, see @ref newExpLeaf, @ref newExpOp and
 * @ref newExpTree.
 */
typedef struct ExpTree {
  /// A field to store a node name or value, as a char array.
  char *data;
  /// The node type impacts requirements for the subtrees and data.
  ExpType type;
  struct ExpTree *left;  ///< The left child/subtree.
  struct ExpTree *right; ///< The right child/subtree.
} ExpTree;

/**
 * @brief The expression tree leaf node constructor.
 * @details All leaves must specify data. This function internally does
 * string duplication, which simplifies most calls to this function.
 * @pre \p name may **not** be NULL.
 * @post The ownership of the data input remains the caller's.
 *
 * @param[in] type The leaf's node type. Is restricted to leaf node types.
 * @param[in] name The leaf node data.
 * @return ExpTree* A newly heap-allocated leaf node.
 */
ExpTree *newExpLeaf(const ExpType type, const char *const name);

/**
 * @brief The expression tree internal (**operator**) node constructor.
 * @details Internal nodes are often operator nodes, meaning they do not
 * need to specify data. Their type is sufficiently informative.
 *
 * The term "operator" refers to mathematical operators, such as addition (+).
 * @pre All pointer arguments must be heap allocated or NULL. Though which
 * subtrees are allowed to be NULL is based on the operator's type.
 * @see ExpType
 * @post Transfers ownership of all heap-allocated arguments to the
 * newly created instance.
 * @post The node data is always NULL.
 *
 * @param[in] type  The node type to assign. Is restricted to internal types.
 * @param[in] left  The left subtree to assign.
 * @param[in] right The right subtree to assign.
 * @return ExpTree* A newly heap-allocated internal node.
 */
ExpTree *newExpOp(const ExpType type, ExpTree *left, ExpTree *right);

/**
 * @brief The expression tree internal (**general**) node constructor.
 * @details This constructor allows all node members to be specified.
 *
 * The term "general" implies that it is possible to construct operator as
 * well as other internal nodes using this constructor. Though, for operator
 * nodes, the @ref newExpOp constructor should be preferred instead.
 * @pre All pointer arguments must be heap allocated or NULL. Though which
 * subtrees are allowed to be NULL is based on the node's type.
 * @see ExpType
 * @post Transfers ownership of all heap-allocated arguments to the
 * newly created instance.
 *
 * @param[in] type  The node type to assign. Is restricted to internal types.
 * @param[in] name  The internal node data.
 * @param[in] left  The left subtree to assign.
 * @param[in] right The right subtree to assign.
 * @return ExpTree* A newly heap-allocated internal node.
 */
ExpTree *newExpTree(const ExpType type, char *name, ExpTree *left,
                    ExpTree *right);

/**
 * @brief Deallocate the given tree recursively.
 * @pre The given tree must not be NULL.
 */
void delExpTree(ExpTree *tree);

/**
 * @brief Print a representation of the given tree to the specified stream.
 * @pre Neither the given \p tree nor stream ( \p where ) may be NULL.
 *
 * @param[in] tree  The tree to print.
 * @param[in] where The stream (destination) to print to.
 */
void printExpTree(const ExpTree *tree, FILE *where);

/**
 * @brief Compute the partial derivative expression w.r.t. the given variable.
 * @details computes \f$ \frac{ \partial f }{ \partial x } \f$ where f
 * is the differentiand and x the differentiation variable.
 *
 * @param[in] expr The expression to derive.
 * @param[in] var  The differentiation variable.
 * @return ExpTree* A newly heap-allocated, partial derivative expression
 * of the input.
 */
ExpTree *derivative(const ExpTree *expr, const char *var);

/**
 * @brief Compute the indefinite integral expression w.r.t. the given variable.
 * @details computes \f$ \int f(x) dx \f$ where f(x) is the integrand and x the
 * integration variable.
 *
 * @param[in] expr The expression to integrate; the integrand.
 * @param[in] var  The integration variable.
 * @return ExpTree* A newly heap-allocated, integral expression of the input.
 */
ExpTree *integral(const ExpTree *expr, const char *var);

/**
 * @brief Compute the definite integral expression w.r.t. the given variable
 * and domain.
 * @details computes \f$ \int_a^b f(x) dx \f$ where f(x) is the integrand and x
 * the integration variable.
 *
 * @param[in] expr   The expression to integrate; the integrand.
 * @param[in] var    The integration variable.
 * @param[in] lowerBound The lower bound of the integration domain.
 * @param[in] upperBound The upper bound of the integration domain.
 * @return ExpTree* A newly heap-allocated, definite integral expression
 * of the input.
 */
ExpTree *definiteIntegral(const ExpTree *expr, const char *var,
                          const ExpTree *lowerBound, const ExpTree *upperBound);

/**
 * @brief Make an exact, deep/recursive copy of the entire expression tree.
 * @details All node data is copied as well.
 *
 * @param[in] src The tree to copy.
 * @return ExpTree* A newly heap-allocated, exact copy.
 */
ExpTree *cpyExpTree(const ExpTree *const src);

/**
 * @brief Verify if the expression is linear.
 *
 * @return true  iff. the expression is linear.
 * @return false else.
 */
bool isLinear(const ExpTree *expr);

/**
 * @brief Verify the exact equality of the given trees.
 * @details Both structure and content (data) must match exactly.
 *
 * Algebraic equivalence is not equality: (1 + x) =\\= (x + 1).
 *
 * Content refers to the data of nodes: (1 + y) =\\= (1 + x).
 *
 * @return true  iff. both trees match exactly.
 * @return false else.
 */
bool isEqual(const ExpTree *expr1, const ExpTree *expr2);

/**
 * @brief Compute the degree of the given monomial expression.
 * @details A monomial is a polynomial of only a single term.
 * The expression can only contain multiplication (*) unary negative (-)
 * and exponentiation (^) with non-negative powers.
 * @pre The expression must be of monomial form.
 *
 * @param[in] expr The monomial expression.
 * @return unsigned int The degree of the expression.
 */
unsigned int degreeMonomial(const ExpTree *expr);

/**
 * @brief A forward declaration of @ref substitute in @ref transformations.h.
 */
ExpTree *substitute(const ExpTree *source, const char *var,
                    const ExpTree *target);

#endif
