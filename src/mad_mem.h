#ifndef MAD_MEM_H
#define MAD_MEM_H

void* mymalloc(char* caller, size_t size);
void* mycalloc(char* caller, size_t nelem, size_t size);
void  myfree(char* rout_name, void* p);

void  mad_mem_handler(int sig);

#endif // MAD_MEM_H

