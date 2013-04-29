#include <stdio.h>
#include <stdlib.h>

#include "mad_err.h"
#include "mad_mem.h"

// public interface

void
mad_mem_handler(int sig)
{
  (void)sig;
  
  puts("+++ memory access outside program range, fatal +++");
  exit(EXIT_FAILURE);
}

