#include <stdio.h>
#include <stdlib.h>

#include "mad_wrap_f.h"
#include "mad_mem.h"

// private constant

#define FREECODE 380226     /* check-code to avoid multiple "free" */
#define MTABLE_SIZE 1000000

// private globals

static char* myfree_caller = "none";

#ifdef _MEM_LEAKS
static int item_no=-1;
static int* mtable[MTABLE_SIZE];
#endif

// public interface

void
mad_mem_handler(int sig)
{
  if (strcmp(myfree_caller, "none") == 0)
    puts("+++ memory access outside program range, fatal +++");
  else
    printf("+++ illegal call to free memory from routine: %s +++\n", myfree_caller);

  puts("good bye");
  exit(EXIT_FAILURE);
}

void*
mymalloc(char* caller, size_t size)
#ifdef _MEM_LEAKS
{
  /* calls malloc, checks for memory granted */
  void* p;
  int* i_p;
  size_t l_size = size + sizeof(double);
  if ((p = malloc(l_size)) == NULL)
  {
    fatal_error("memory overflow, called from routine:", caller);
  } 
  i_p = (int*) p;
  mtable[-item_no]=i_p;
  *i_p++ = item_no;
  *i_p = l_size;
  fprintf(stderr,"ALLOCATE called by %s \n",caller);
  fprintf(stderr,"[Allocated item %i (size %i)] \n",item_no,l_size);
  if (item_no == -MTABLE_SIZE)
  {
    fatal_error("Too many allocs!!!", "MTABLE_SIZE");
  }
  item_no = item_no - 1;
  return (void *)((char*)p+sizeof(double));
}
#endif
#ifndef _MEM_LEAKS
{
  /* calls malloc, checks for memory granted */
  void* p;
  int* i_p;
  size_t l_size = size + sizeof(double)+2;
/*  printf("xxxx %d xxxx\n",l_size);*/
  if ((p = malloc(l_size)) == NULL)
  {
    fatal_error("memory overflow, called from routine:", caller);
  } 
  i_p = (int*) p; *i_p = FREECODE;
  return (void *)((char*)p+sizeof(double));
}
#endif

void*
mycalloc(char* caller, size_t nelem, size_t size)
#ifdef _MEM_LEAKS
{
  /* calls calloc, checks for memory granted */
  void* p;
  int* i_p;
  size_t l_size = nelem*size + sizeof(double);
  if ((p = calloc(1, l_size)) == NULL)
  {
    fatal_error("memory overflow, called from routine:", caller);
  } 
  mtable[-item_no]=p;
  i_p = (int*) p;
  *i_p++ = item_no;
  *i_p = l_size;
  fprintf(stderr,"ALLOCATE called by %s \n",caller);
  fprintf(stderr,"[Allocated item %i (size %i)] \n",item_no,l_size);
  if (item_no == -MTABLE_SIZE)
  {
    fatal_error("Too many allocs!!!", "MTABLE_SIZE");
  }
  item_no = item_no - 1;
  return (void *)((char*)p+sizeof(double));
}
#endif
#ifndef _MEM_LEAKS
{
  /* calls calloc, checks for memory granted */
  void* p;
  int* i_p;
  size_t l_size = nelem*size + sizeof(double);
  if ((p = calloc(1, l_size)) == NULL)
  {
    fatal_error("memory overflow, called from routine:", caller);
  } 
  i_p = (int*) p; *i_p = FREECODE;
  return ((char*)p+sizeof(double));
}
#endif

void
myfree(char* rout_name, void* p)
#ifdef _MEM_LEAKS
{
  int my_size,my_item_no,myend,old_item_no;
  char* l_p = (char*)p - sizeof(double);
  int* i_p = (int*) l_p;
  myfree_caller = rout_name;
/* Look for the integer address (backwards) in mtable */
  myend=-item_no-1;
  old_item_no=0;
  while (myend > 0)
  {
    if (mtable[myend] == i_p)
    {
      old_item_no=-myend;
      mtable[myend]=NULL;
      break;
    }
    else
    {
      myend=myend-1;
    }
  }
  if ( old_item_no == 0)
  {
    fatal_error("Free memory error!!!, called from routine:", myfree_caller);
  }
  my_item_no = *i_p++;
  my_size = *i_p;
  if (my_item_no != old_item_no)
  {
    fatal_error("Memory item number discrepancy!!!:", myfree_caller);
  }
  fprintf(stderr,"DEALLOCATE called by %s \n",myfree_caller);
  fprintf(stderr,"[Deallocated item %i (size %i)] \n",my_item_no,my_size);
  free(l_p);
  myfree_caller = "none";
}
#endif
#ifndef _MEM_LEAKS
{
  char* l_p = (char*)p - sizeof(double);
  int* i_p = (int*) l_p;
  myfree_caller = rout_name;
  if (*i_p == FREECODE)
  {
    *i_p = 0; free(l_p);
  }
  myfree_caller = "none";
}
#endif

