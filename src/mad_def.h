#ifndef MAD_DEF_H
#define MAD_DEF_H

// Temporary file:
// these constants should be split over their respective modules...
// and most of them should be replaced by constant variables or enums

/* preparation of Touschek */
/* defined constants for word lengths etc. */
#define ALIGN_MAX 14        /* alignment error array length */
#define EFIELD_TAB 42       /* field error array length for ESAVE table */
#define FIELD_MAX 42        /* field error array length */
#define RFPHASE_MAX 42      /* rf-phase error array length */
#define SEQ_DUMP_LEVEL 0    /* chooses amount of dumped output */
#define NAME_L 48           /* internal name length */
#define TITLE_SIZE 114      /* Size of the title for gnuplot ploting in tracking mode */
#define PTC_NAMES_L 13      /* Number of ptc variables treated in select_ptc_normal */
#define MAX_ROWS 101        /* Initial size of ptc_normal table */
#define FNAME_L 240         /* for file names */
#define INVALID 1.e20       /* used for erroneous value requests */
#define AUX_LG 50000        /* initial size for ancillary buffers */
#define MAX_ITEM  1000      /* initial # of items in tok_list etc. */
#define MAX_D_ITEM 30000    /* initial storage size for doubles */
#define MAX_LINE 20000      /* max. input line length (has to stay fixed) */
#define MAX_LOOP 100        /* max. count for (possibly circular) calls */
#define MAX_COND 100        /* max. nesting level for "if" and "while" */
#define MAX_TYPE 11         /* for SXF output */
#define MAX_TAG 50          /* for SXF output */
#define CHAR_BUFF_SIZE 100000 /* size of dynamic char_buff members */
#define IN_BUFF_SIZE 500000 /* initial size of buffer for command groups */
#define LINE_FILL 240        /* max. line length -2 for "save" output */
#define LINE_F_MAD8 70      /* the same, for mad-8 format */
#define MADX_LINE_MAX 78         /* for SXF output */
#define MATCH_WORK 10       /* no. of work spaces in matching */
#define USER_TABLE_LENGTH 100 /* initial length of user defined tables */
#define MAXARRAY 1000       /* max. length of apex tables in aperture module*/
#define DQ_DELTAP 1.e-6     /* deltap for difference calculation of chrom. */

#define E_D_MAX 500         /* max. length of extra displacement tables (per element) */
#define E_D_LIST_CHUNK 1000  /* chunk to allocate memory for extra displacement tables */

#define MADX_LONG      1
#define MADX_DOUBLE    2
#define MADX_STRING    3

#define MAX_TFS_ROW 2000  /* max. number of rows for SDDS  conversion */
#define MAX_TFS_COL 500   /* max. number of columns for SDDS conversion */

#endif // MAD_DEF_H

