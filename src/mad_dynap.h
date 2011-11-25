#ifndef MAD_DYNAP_H
#define MAD_DYNAP_H

// types

struct in_cmd;

// interface

void dynap_tables_create(struct in_cmd* cmd);
void track_dynap(struct in_cmd* cmd);

#endif // MAD_DYNAP_H

