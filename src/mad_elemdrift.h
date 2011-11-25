#ifndef MAD_ELEMDRIFT_H
#define MAD_ELEMDRIFT_H

// types

struct node;

// interface

struct element* get_drift(double length);
int   add_drifts(struct node* c_node, struct node* end);

#endif // MAD_ELEMDRIFT_H
