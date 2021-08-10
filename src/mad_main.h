#ifndef MAD_MAIN_H
#define MAD_MAIN_H

// readonly information about program's command line arguments and stack base
extern int                const mad_argc;
extern const char* const* const mad_argv;
extern const void* const        mad_stck;
extern char rel_path_dir[100];

// public interface to run MADX as a library
void mad_init(int argc, char *argv[]);
void mad_fini(void);
void mad_run (void);

#endif // MAD_MAIN_H

