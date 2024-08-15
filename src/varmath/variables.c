#include "variables.h"

Domain *newDomain(char *var, Interval domain) {
  assert(var != NULL);

  Domain *list = (Domain *)malloc(sizeof(Domain));
  list->var = var;
  list->domain = domain;
  list->next = NULL;
  return list;
}

Domain *appDomainElem(Domain *tail, Domain *head) {
  assert(head != NULL);
  assert(head->next == NULL);

  head->next = tail;
  return head;
}

Domain *newDomainElem(Domain *tail, char *var, Interval domain) {
  Domain *head = newDomain(var, domain);
  return appDomainElem(tail, head);
}

void delDomain(Domain *list) {
  if (list->next != NULL)
    delDomain(list->next);
  assert(list->var != NULL);
  free(list->var);
  free(list);
}

void printDomain(Domain *list, FILE *where) {
  assert(list->var != NULL);
  fprintf(where, "%s in ", list->var);

  printInterval(&list->domain, where);

  fprintf(where, "; ");
  if (list->next != NULL)
    printDomain(list->next, where);
}