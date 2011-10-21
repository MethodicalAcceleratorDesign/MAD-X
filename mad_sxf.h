#ifndef MAD_SXF_H
#define MAD_SXF_H

void  get_sxf_names(void);
void  pro_sxf(struct in_cmd* cmd);
void  sxf_read(struct command* comm);
int   sxf_decin(char* p, int count); /* decode one SXF input item, store */
void  sxf_fill_command(struct command* comm, int ntok, char** toks);
int   sxf_align_fill(int start, int end, int ntok, char** toks, double* vec);
int   sxf_field_fill(int start, int end, int ntok, char** toks, double* vec);
void  sxf_body_fill(struct command* comm, int start, int end, int ntok, char** toks, double length);
void  sxf_write(struct command* comm, FILE* out);
void  sxf_rtag(void);

#endif // MAD_SXF_H
