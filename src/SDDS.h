 /*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: SDDS.h
 * purpose: define structures and prototypes for use with SDDS
 *          format
 *
 * Michael Borland, 1993
 $Log: not supported by cvs2svn $
 Revision 1.85  2009/06/02 18:17:15  soliday
 Updated so that the MPI functions will be visable from shared libraries.

 Revision 1.84  2009/06/01 20:57:55  soliday
 Updated to work with SDDS_MPI_IO=0

 Revision 1.83  2009/05/08 16:05:18  shang
 added MPI non native writing routines

 Revision 1.82  2008/11/11 21:05:18  soliday
 Updated to fix an issue on vxWorks.

 Revision 1.81  2008/09/18 20:15:07  soliday
 Added SDDS_SetColumnFromFloats and SDDS_GetColumnInFloats

 Revision 1.80  2008/08/04 14:04:57  shang
 added collective_io to MPI_DATASET structure

 Revision 1.79  2008/07/21 14:45:28  shang

 redefined SDDS_MPI_TotalRowsCount

 Revision 1.78  2008/04/14 18:36:57  shang
 now includes the definition of parallel SDDS.

 Revision 1.77  2007/10/12 18:03:44  shang
 added titleBuffer member to SDDS_DATASET

 Revision 1.76  2007/09/27 14:25:53  shang
 added SDDS_Parallel_InitializeOutput

 Revision 1.75  2006/12/07 16:40:54  soliday
 Added SDDS_GetToken2 and SDDS_ScanData2 which are now used when reading
 ascii array data.

 Revision 1.74  2006/09/15 18:15:48  borland
 Added prototype for SDDS_CheckEndOfFile.

 Revision 1.73  2006/08/31 15:06:55  soliday
 Updated to work with SDDS2

 Revision 1.72  2006/05/25 21:23:55  shang
 added SDDS_GetColumnInShort

 Revision 1.71  2006/03/16 18:50:02  soliday
 Added the SDDS_BreakIntoLockedFile function.

 Revision 1.70  2005/11/04 22:46:59  soliday
 Updated code to be compiled by a 64 bit processor.

 Revision 1.69  2005/01/13 16:52:51  shang
 added is_string argument to is_memory() function to get the memory data type.

 Revision 1.68  2004/11/15 16:59:19  soliday
 Moved and renamed the procedures to write non native binary data from
 sddsendian to SDDSlib.

 Revision 1.67  2004/05/13 18:16:08  soliday
 Added SDDS_AppendToArrayVararg

 Revision 1.66  2004/02/27 21:57:57  soliday
 Modified SDDS_SetRowsOfInterest

 Revision 1.65  2004/02/27 21:14:08  soliday
 Modified SDDS_RegisterProgramName

 Revision 1.64  2004/02/05 22:56:15  soliday
 Added missing declaration for SDDS_CheckTabularData

 Revision 1.63  2003/10/17 16:11:29  soliday
 Added SDDS_VerifyColumnExists, SDDS_VerifyParameterExists,
 and SDDS_VerifyArrayExists

 Revision 1.62  2003/09/23 21:49:23  soliday
 Added the SDDS_MATCH_EXCLUDE_STRING definition.

 Revision 1.61  2003/09/11 22:31:35  soliday
 Added SDDS_GetParameterAsLong and SDDS_ConvertToLong

 Revision 1.60  2003/07/22 20:01:16  soliday
 Added support for Kylix.

 Revision 1.59  2003/02/15 22:54:16  borland
 Added prototype for SDDS_DoFSync().

 Revision 1.58  2002/08/14 15:40:12  soliday
 Added Open License

 Revision 1.57  2002/07/25 00:10:37  borland
 Added definition of SDDS_PRINT_BUFLEN.

 Revision 1.56  2002/07/24 20:23:13  shang
 added SDDS_InitializeInputFromSearchPath()

 Revision 1.55  2002/06/15 02:28:00  borland
 Added SDDS_ForceActive().

 Revision 1.54  2002/06/05 16:20:41  shang
 added SDDS_GetParameterAsString()

 Revision 1.53  2002/05/26 02:07:12  borland
 Added autoRecovered variable to SDDS_DATASET to allow recording fact that
 auto-read-recovery has occurred.

 Revision 1.52  2002/02/18 19:27:47  soliday
 Added the ability to enable or disable fsync.

 Revision 1.51  2002/01/15 22:46:05  soliday
 Added SetRowCountMode, UpdateRowCount, SetAutoReadRecovery

 Revision 1.50  2002/01/11 17:08:19  soliday
 Added SDDS_SyncDataSet

 Revision 1.49  2002/01/07 21:33:16  borland
 Added prototypes for SDDS_GetRowLimit() and SDDS_SetRowLimit().

 Revision 1.48  2002/01/07 19:31:23  soliday
 Fixed declaration for SDDS_SetNoRowCounts.

 Revision 1.47  2001/12/22 21:21:32  soliday
 Fixed RW_ASSOCIATES typo.

 Revision 1.46  2001/12/21 03:29:39  soliday
 Fixed some declarations.

 Revision 1.45  2001/11/30 15:35:13  borland
 Changes by H. Shang: addition of SDDS_GotoPage() and changes required to
 support this.

 Revision 1.44  2001/02/19 18:01:25  borland
 Added prototype for SDDS_GetParameterByIndex().

 Revision 1.43  2000/11/27 17:15:24  borland
 Fixed typo in previous change.

 Revision 1.42  2000/11/27 17:06:27  borland
 Added prototype for SDDS_FreeStringData().

 Revision 1.41  2000/10/31 21:24:53  soliday
 Fixed the declarations of some new functions so that they work with WIN32.

 Revision 1.40  2000/10/24 01:16:45  borland
 Added prototype for SDDS_GetArrayInDoubles().

 Revision 1.39  2000/09/19 20:14:11  borland
 Added function prototypes and changed some structures (SDDS_DATASET and
 SDDS_LAYOUT) to support reading of little/big endian files.

 Revision 1.38  2000/06/20 18:02:42  borland
 Added prototype for SDDS_SetDefaultIOBufferSize.

 Revision 1.37  2000/06/20 14:28:37  soliday
 Moved some definitions from SDDS_internal.h to SDDS.h so that sddsendian
 can access them.

 Revision 1.36  2000/04/14 16:36:12  soliday
 gzFile now is defined even when we are not using the zLib library.

 Revision 1.35  2000/04/13 17:04:58  soliday
 Fixed prototype for SDDS_NumberOfErrors and SDDS_ClearErrors

 Revision 1.34  2000/03/30 19:58:35  borland
 Added SDDS_DisconnectFile() and SDDS_ReconnectFile(), plus disconnected member
 to SDDS_LAYOUT structure.

 Revision 1.33  2000/03/30 17:26:57  soliday
 Added SDDS_InitializeAppendToPage

 Revision 1.32  2000/03/27 20:25:47  borland
 Added prototype for SDDS_ShortenTable().

 Revision 1.31  2000/02/25 23:15:12  evans
 Added SDDS_Free.  Memory allocated with SDDS_Malloc and others should
 be freed with SDDS_Free.  It is necessary for WIN32 when degugging
 SDDS applications using the release version of the SDDS library
 because the memory allocation models for debug and release are
 inconsistent.

 Revision 1.30  2000/01/18 20:06:46  soliday
 Added support for ZLIB.

 Revision 1.29  1999/11/03 22:42:47  soliday
 Added WriteBinaryString declaration

 Revision 1.28  1999/09/16 22:03:15  soliday
 Modified the way global variables are exported and imported to a WIN32 dll

 Revision 1.27  1999/09/14 18:03:28  soliday
 Added export commands for WIN32 dll files.

 Revision 1.26  1999/06/01 19:17:05  borland
 Added SDDS_Calloc.

 Revision 1.25  1999/04/14 13:58:47  borland
 Fixed some function types and returns.

 Revision 1.24  1999/02/04 21:13:19  soliday
 Added prototype for SDDS_GetColumnsInLong

 Revision 1.23  1998/12/16 21:20:35  borland
 Modified prototypes for SDDS_TransferAll* functions and added required
 macros.

 Revision 1.22  1998/08/25 01:45:47  borland
 Moved SDDS_FILEBUFFER_SIZE definition here from SDDS_binary.c

 Revision 1.21  1998/06/30 21:30:48  borland
 Added prototype for SDDS_FreeStringArray().

 Revision 1.20  1997/12/19 16:55:45  borland
 Fixed SDDS_RowCount macro (more parentheses).  Added prototype for
 SDDS_Malloc.  Added new "type": SDDS_ANY_FLOATING_TYPE.

 Revision 1.19  1997/07/08 14:41:17  borland
 Added SDDS_RowCount macro.

 * Revision 1.18  1996/07/05  16:36:02  borland
 * Modified prototypes for SDDS_PrintTypedValue and SDDS_SprintTypedValue
 * to include new mode argument.
 *
 * Revision 1.17  1996/03/27  17:51:08  borland
 * Added prototype and macros for SDDS_GetErrorMessages().
 *
 * Revision 1.16  1996/03/21  20:07:32  borland
 * Added mode argument to prototype for SDDS_ReadPageSparse().
 *
 * Revision 1.15  1996/03/21  16:39:50  borland
 * Added definitions for major revision of SDDS library.  Added file buffer
 * structure to SDDS_DATASET.  Added prototype for SDDS_ReadPageSparse().
 *
 * Revision 1.14  1996/02/12  17:18:55  borland
 * Added prototype for SDDS_SetAutoCheckMode() and mode flags to go with it.
 *
 * Revision 1.13  1996/01/20  05:38:03  borland
 * Added prototypes for SDDS_AppendLayout().
 *
 * Revision 1.12  1996/01/19  00:18:09  borland
 * SDDS.h: Added popenUsed to SDDS_LAYOUT structure.
 * scan.h: Added prototypes for unpacking functions.
 *
 * Revision 1.11  1996/01/18  04:31:56  borland
 * Added prototype for SDDS_GetInternalColumn()
 *
 * Revision 1.10  1996/01/14  01:56:27  borland
 * Added prototype and flag definitions for SDDS_FilterByNumScan().
 *
 * Revision 1.9  1996/01/10  17:23:57  borland
 * Added prototype for SDDS_LockFile().
 *
 * Revision 1.8  1996/01/09  23:22:58  borland
 * Added prototypes for SDDS_TransferAll*Definitions routines.
 *
 * Revision 1.7  1995/12/12  10:02:28  borland
 * Added prototype for SDDS_SetNoRowCounts().
 *
 * Revision 1.6  1995/12/12  04:30:38  borland
 * Added layout_written field to SDDS_LAYOUT struct; added  file_had_data
 * field to SDDS_DATASET struct; both in support of write/append to
 * files with no_row_counts=1.
 *
 * Revision 1.5  1995/11/13  16:08:13  borland
 * Added first_row_in_mem and writing_page entries to SDDS_DATASET structure
 * to support partially flushed tabular data.
 *
 * Revision 1.4  1995/11/08  04:25:53  borland
 * Added prototypes for new routines FreeArray, DefineParameterLikeArray,
 * and DefineColumnLikeArray.
 *
 * Revision 1.3  1995/09/12  03:18:41  borland
 * Added flag bit #define's for name validity control.
 * Added SDDS_CI_* #defines for SDDS_SetRowsOfInterest case-insensitive
 *   mode.
 * Added SDDS_NOCASE_COMPARE for SDDS_MatchRowsOfInterest case-insenstive
 *   mode.
 *
 * Revision 1.2  1995/09/05  21:14:46  saunders
 * First test release of the SDDS1.5 package.
 *
 */

#if !defined(_SDDS_)
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#define SDDS_READMODE 1
#define SDDS_WRITEMODE 2
#define SDDS_MEMMODE 3

#define SDDS_MPI_READ_ONLY 0x0001UL
#define SDDS_MPI_WRITE_ONLY 0x0002UL
#define SDDS_MPI_READ_WRITE  0x0004UL
#define SDDS_MPI_STRING_COLUMN_LEN 16

#if defined(SDDS_MPI_IO) && SDDS_MPI_IO
#include "mpi.h"
#include "mdb.h"

extern char SDDS_mpi_error_str[MPI_MAX_ERROR_STRING];
extern int32_t SDDS_mpi_error_str_len;
extern char *SDDS_MPI_FILE_TYPE[];

typedef struct {
  MPI_File	        MPI_file;	/*MPIO file handle			*/
  MPI_Comm	        comm;		/*communicator				*/
  MPI_Info	        File_info;	/*file information			*/
  int32_t               myid;           /* This process's rank                  */
  int32_t               n_processors;   /* Total number of processes        */
  MPI_Offset            file_offset, file_size, column_offset;  /* number of bytes in one row for row major order*/
  short                 collective_io;
  int32_t               n_page;         /* index of current page*/
  int32_t               n_rows;         /* number of rows that current processor holds */
  int32_t               total_rows;     /* the total number of rows that all processor hold*/
  int32_t               end_of_file;    /* flag for end of MPI_file */
  int32_t               master_read;   /*determine if master processor read the page data or not*/
  int32_t               start_row, end_row; /* the start row and end row that current processor's data that is going to be written to output or read from input */
} MPI_DATASET;
#endif

#if defined(zLib)
#include "zlib.h"
#endif

#ifdef __cplusplus 
extern "C" {
#endif


#define _SDDS_ 1

#define SDDS_VERSION 2

#include "SDDStypes.h"
#if defined(_WIN32)
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#define PRId32 "ld"
#define SCNd32 "ld"
#define PRIu32 "lu"
#define SCNu32 "lu"
#define INT32_MAX (2147483647)
#else
#if defined(vxWorks)
#define PRId32 "ld"
#define SCNd32 "ld"
#define PRIu32 "lu"
#define SCNu32 "lu"
#define INT32_MAX (2147483647)
#else
#include <inttypes.h>
#endif
#endif

#undef epicsShareFuncSDDS
#undef epicsShareExtern
#if (defined(_WIN32) && !defined(__CYGWIN32__)) || (defined(__BORLANDC__) && defined(__linux__))
#if defined(EXPORT_SDDS)
#define epicsShareFuncSDDS  __declspec(dllexport)
#define epicsShareExtern extern __declspec(dllexport)
#else
#define epicsShareFuncSDDS
#define epicsShareExtern extern __declspec(dllimport)
#endif
#else
#define epicsShareFuncSDDS
#define epicsShareExtern extern
#endif


#define RW_ASSOCIATES 0

#define SDDS_BINARY 1
#define SDDS_ASCII  2
#define SDDS_PARALLEL 3
#define SDDS_NUM_DATA_MODES 2

  /*
extern char *SDDS_type_name[SDDS_NUM_TYPES];
extern int32_t SDDS_type_size[SDDS_NUM_TYPES];
extern char *SDDS_data_mode[SDDS_NUM_DATA_MODES];
  */

epicsShareExtern char *SDDS_type_name[SDDS_NUM_TYPES];
epicsShareExtern int32_t SDDS_type_size[SDDS_NUM_TYPES];
extern char *SDDS_data_mode[SDDS_NUM_DATA_MODES];

/* this shouldn't be changed without changing buffer sizes in namelist routines: */
#define SDDS_MAXLINE 1024

#define SDDS_PRINT_BUFLEN 16834

#define SDDS_NORMAL_DEFINITION 0
#define SDDS_WRITEONLY_DEFINITION 1

typedef struct {
    char *name, *symbol, *units, *description, *format_string;
    int32_t type;
    char *fixed_value;
    /* these are used internally and are not set by the user: */
    int32_t definition_mode, memory_number;
    } PARAMETER_DEFINITION;
#define SDDS_PARAMETER_FIELDS 7

typedef struct {
    char *name, *symbol, *units, *description, *format_string;
    int32_t type, field_length;
    /* these are used internally and are not set by the user: */
    int32_t definition_mode, memory_number, pointer_number;
    } COLUMN_DEFINITION;
#define SDDS_COLUMN_FIELDS 7

typedef struct {
    char *name, *symbol, *units, *description, *format_string, *group_name;
    int32_t type, field_length, dimensions;
    } ARRAY_DEFINITION;
#define SDDS_ARRAY_FIELDS 9

typedef struct {
    char *name, *filename, *path, *description, *contents;
    int32_t sdds;
    } ASSOCIATE_DEFINITION;
#define SDDS_ASSOCIATE_FIELDS 6

typedef struct {
  int32_t index;
  char *name;
} SORTED_INDEX;
int SDDS_CompareIndexedNames(void *s1, void *s2);
int SDDS_CompareIndexedNamesPtr(const void *s1, const void *s2);

typedef struct {
    int32_t mode, lines_per_row, no_row_counts, fixed_row_count, fsync_data;
    int32_t additional_header_lines;
    } DATA_MODE;

#if !defined(zLib)
  typedef void *voidp;
  typedef voidp gzFile;
#endif

typedef struct {
    int32_t n_columns, n_parameters, n_associates, n_arrays;

    char *description;
    char *contents;
    int32_t version;
    short layout_written;

    DATA_MODE data_mode;
    COLUMN_DEFINITION *column_definition;
    PARAMETER_DEFINITION *parameter_definition;
    ARRAY_DEFINITION *array_definition;
    ASSOCIATE_DEFINITION *associate_definition;
    SORTED_INDEX **column_index, **parameter_index, **array_index;

    char *filename;
    FILE *fp;
    gzFile *gzfp;
    short disconnected;
    short gzipFile;
    short popenUsed;
    uint32_t byteOrderDeclared;
    } SDDS_LAYOUT;

typedef struct {
    ARRAY_DEFINITION *definition;
    int32_t *dimension, elements;
    /* Address of array of data values, stored contiguously.
     * For STRING data, the "data values" are actually the addresses of the individual strings.
     */
    void *data;
    void *pointer;
    } SDDS_ARRAY;

typedef struct {
  char *data, *buffer;
  int32_t bytesLeft, bufferSize;
} SDDS_FILEBUFFER ;

#define SDDS_FILEBUFFER_SIZE  262144

typedef struct {
  SDDS_LAYOUT layout, original_layout;
  short swapByteOrder;
  SDDS_FILEBUFFER fBuffer;
  SDDS_FILEBUFFER titleBuffer;
    int32_t page_number;

    short mode; /*file mode*/
    short page_started; /*0 or 1 for page not started or already started*/
    int32_t pages_read; /*the number of pages read so far*/
    int32_t endOfFile_offset;/*the offset in the end of the file*/
    int32_t *pagecount_offset; /*the offset of each read page */ 
    int32_t rowcount_offset;  /* ftell() position of row count */
    int32_t n_rows_written;   /* number of tabular data rows written to disk */
    int32_t last_row_written; /* virtual index of last row written */
    int32_t first_row_in_mem; /* virtual index of first row in memory */
    short writing_page;    /* 1/0 for writing/not-writing page */
    int32_t n_rows_allocated; /* number of tabular data rows for which space is alloc'd */
    int32_t n_rows;           /* number of rows stored in memory */
    int32_t *row_flag;        /* accept/reject flags for rows stored in memory */
    short file_had_data;   /* indicates that file being appended to already had some data in it (i.e.,
                            * more than just a header.  Affects no_row_counts=1 output.
                            */
    short autoRecover;
    short autoRecovered;
  short parallel_io;        /*flag for parallel SDDS */
    int32_t n_of_interest;
    int32_t *column_order;          /* column_order[i] = internal index of user's ith column */
    int32_t *column_flag;           /* column_flag[i] indicates whether internal ith column has been selected */

    /* array of SDDS_ARRAY structures for storing array data */
    SDDS_ARRAY *array;

    /* array for parameter data.  The address of the data for the ith parameter
     * is parameter[i].  So *(<type-name> *)parameter[i] gives the data itself.  For type
     * SDDS_STRING the "data itself" is actually the address of the string, and the type-name
     * is "char *".
     */
    void **parameter;

    /* array for tabular data.  The starting address for the data for the ith column
     * is data[i].  So ((<type-name> *)data[i])[j] gives the jth data value.  
     * For type SDDS_STRING the "data value" is actually the address of the string, and
     * the type-name is "char *".
     */
    void **data;
  short column_major;
#if defined(SDDS_MPI_IO) && SDDS_MPI_IO
  MPI_DATASET *MPI_dataset;
#endif
    } SDDS_DATASET;

typedef SDDS_DATASET SDDS_TABLE;

/* prototypes for routines to prepare and write SDDS files */
epicsShareFuncSDDS extern int32_t SDDS_InitializeOutput(SDDS_DATASET *SDDS_dataset, int32_t data_mode,
						     int32_t lines_per_row, char *description,
						     char *contents, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_Parallel_InitializeOutput(SDDS_DATASET *SDDS_dataset, char *description,
						     char *contents, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_InitializeAppend(SDDS_DATASET *SDDS_dataset, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_InitializeAppendToPage(SDDS_DATASET *SDDS_dataset, char *filename, 
                                                           int32_t updateInterval,
							   int32_t *rowsPresentReturn);
epicsShareFuncSDDS extern int32_t SDDS_DisconnectFile(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ReconnectFile(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_FreeStringData(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_Terminate(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern void SDDS_SetTerminateMode(uint32_t mode);
#define TERMINATE_DONT_FREE_TABLE_STRINGS 0x0001
#define TERMINATE_DONT_FREE_ARRAY_STRINGS 0x0002
epicsShareFuncSDDS extern int32_t SDDS_SetRowCountMode(SDDS_DATASET *SDDS_dataset, uint32_t mode);
#define SDDS_VARIABLEROWCOUNT 0x0001UL
#define SDDS_FIXEDROWCOUNT 0x0002UL
#define SDDS_NOROWCOUNT 0x0004UL
epicsShareFuncSDDS extern int32_t SDDS_SetAutoReadRecovery(SDDS_DATASET *SDDS_dataset, uint32_t mode);
#define SDDS_NOAUTOREADRECOVER 0x0001UL
#define SDDS_AUTOREADRECOVER 0x0002UL
epicsShareFuncSDDS extern int32_t SDDS_UpdateRowCount(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern void SDDS_DisableFSync(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern void SDDS_EnableFSync(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_DoFSync(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_DefineParameter(SDDS_DATASET *SDDS_dataset, char *name, char *symbol, char *units, char *description, 
                           char *format_string, int32_t type, char *fixed_value);
epicsShareFuncSDDS extern int32_t SDDS_DefineParameter1(SDDS_DATASET *SDDS_dataset, char *name, char *symbol, char *units, char *description, 
                           char *format_string, int32_t type, void *fixed_value);
epicsShareFuncSDDS extern int32_t SDDS_DefineColumn(SDDS_DATASET *SDDS_dataset, char *name, char *symbol, char *units, char *description, 
                        char *format_string, int32_t type, int32_t field_length);
// epicsShareFuncSDDS extern int32_t SDDS_DefineArray(SDDS_DATASET *SDDS_dataset, char *name, char *symbol, char *units, char *description, 
epicsShareFuncSDDS extern int32_t SDDS_DefineArray(SDDS_DATASET *SDDS_dataset, const char *name, char *symbol, char *units, char *description, 
                        char *format_string, int32_t type, int32_t field_length, int32_t dimensions, char *group_name);
epicsShareFuncSDDS extern int32_t SDDS_DefineAssociate(SDDS_DATASET *SDDS_dataset, char *name,
                                 char *filename, char *path, char *description, char *contents, int32_t sdds);
epicsShareFuncSDDS extern int32_t SDDS_IsValidName(char *name, char *dataClass);
epicsShareFuncSDDS extern int32_t SDDS_SetNameValidityFlags(uint32_t flags);
#define SDDS_ALLOW_ANY_NAME 0x0001UL
#define SDDS_ALLOW_V15_NAME 0x0002UL
epicsShareFuncSDDS extern int32_t SDDS_DefineSimpleColumn(SDDS_DATASET *SDDS_dataset, char *name, char *unit, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_DefineSimpleParameter(SDDS_DATASET *SDDS_dataset, char *name, char *unit, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_DefineSimpleColumns(SDDS_DATASET *SDDS_dataset, int32_t number, char **name, char **unit, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_DefineSimpleParameters(SDDS_DATASET *SDDS_dataset, int32_t number, char **name, char **unit, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_SetNoRowCounts(SDDS_DATASET *SDDS_dataset, int32_t value);
epicsShareFuncSDDS extern int32_t SDDS_WriteLayout(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_EraseData(SDDS_DATASET *SDDS_dataset);

epicsShareFuncSDDS extern int32_t SDDS_ProcessColumnString(SDDS_DATASET *SDDS_dataset, char *string, int32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_ProcessParameterString(SDDS_DATASET *SDDS_dataset, char *string, int32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_ProcessArrayString(SDDS_DATASET *SDDS_dataset, char *string);
epicsShareFuncSDDS extern int32_t SDDS_ProcessAssociateString(SDDS_DATASET *SDDS_dataset, char *string);

epicsShareFuncSDDS extern int32_t SDDS_InitializeCopy(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source, char *filename, char *filemode);
epicsShareFuncSDDS extern int32_t SDDS_CopyLayout(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_AppendLayout(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source, uint32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_CopyPage(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
#define SDDS_CopyTable(a, b) SDDS_CopyPage(a, b)
epicsShareFuncSDDS extern int32_t SDDS_CopyParameters(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_CopyArrays(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_CopyColumns(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_CopyRowsOfInterest(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);
epicsShareFuncSDDS extern int32_t SDDS_CopyRow(SDDS_DATASET *SDDS_target, int32_t target_row, SDDS_DATASET *SDDS_source, int32_t source_srow);
epicsShareFuncSDDS extern int32_t SDDS_CopyRowDirect(SDDS_DATASET *SDDS_target, int32_t target_row, SDDS_DATASET *SDDS_source, int32_t source_row);
epicsShareFuncSDDS extern int32_t SDDS_CopyAdditionalRows(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source);

epicsShareFuncSDDS extern void SDDS_DeferSavingLayout(int32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_SaveLayout(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_RestoreLayout(SDDS_DATASET *SDDS_dataset);

#define SDDS_BY_INDEX 1
#define SDDS_BY_NAME  2

epicsShareFuncSDDS extern int32_t SDDS_StartPage(SDDS_DATASET *SDDS_dataset, int32_t expected_n_rows);
#define SDDS_StartTable(a, b) SDDS_StartPage(a, b)
epicsShareFuncSDDS extern int32_t SDDS_ClearPage(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_LengthenTable(SDDS_DATASET *SDDS_dataset, int32_t n_additional_rows);
epicsShareFuncSDDS extern int32_t SDDS_ShortenTable(SDDS_DATASET *SDDS_dataset, int32_t rows);
#define SDDS_SET_BY_INDEX SDDS_BY_INDEX
#define SDDS_SET_BY_NAME SDDS_BY_NAME
#define SDDS_PASS_BY_VALUE 4
#define SDDS_PASS_BY_REFERENCE 8
#define SDDS_PASS_BY_STRING 16
epicsShareFuncSDDS extern int32_t SDDS_SetParameters(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetRowValues(SDDS_DATASET *SDDS_dataset, int32_t mode, int32_t row, ...);
epicsShareFuncSDDS extern int32_t SDDS_WritePage(SDDS_DATASET *SDDS_dataset);
#define SDDS_WriteTable(a) SDDS_WritePage(a)
epicsShareFuncSDDS extern int32_t SDDS_UpdatePage(SDDS_DATASET *SDDS_dataset, uint32_t mode);
#define FLUSH_TABLE 0x1UL
#define SDDS_UpdateTable(a) SDDS_UpdatePage(a, 0)
epicsShareFuncSDDS extern int32_t SDDS_SyncDataSet(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_SetColumn(SDDS_DATASET *SDDS_dataset, int32_t mode, void *data, int32_t rows, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetColumnFromDoubles(SDDS_DATASET *SDDS_dataset, int32_t mode, double *data, int32_t rows, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetColumnFromFloats(SDDS_DATASET *SDDS_dataset, int32_t mode, float *data, int32_t rows, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetColumnFromLongs(SDDS_DATASET *SDDS_dataset, int32_t mode, int32_t *data, int32_t rows, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetParametersFromDoubles(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);

#define SDDS_GET_BY_INDEX SDDS_BY_INDEX
#define SDDS_GET_BY_NAME SDDS_BY_NAME
epicsShareFuncSDDS extern int32_t SDDS_GetColumnInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_GetParameterInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_GetArrayInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_GetAssociateInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_ChangeColumnInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_ChangeParameterInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_ChangeArrayInformation(SDDS_DATASET *SDDS_dataset, char *field_name, void *memory, int32_t mode, ...);

epicsShareFuncSDDS extern void SDDS_SetReadRecoveryMode(int32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_SetDefaultIOBufferSize(int32_t bufferSize);

/* prototypes for routines to read and use SDDS files  */
epicsShareFuncSDDS extern int32_t SDDS_InitializeInputFromSearchPath(SDDS_DATASET *SDDSin, char *file);
epicsShareFuncSDDS extern int32_t SDDS_InitializeInput(SDDS_DATASET *SDDS_dataset, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_ReadLayout(SDDS_DATASET *SDDS_dataset, FILE *fp);
epicsShareFuncSDDS extern int32_t SDDS_InitializeHeaderlessInput(SDDS_DATASET *SDDS_dataset, char *filename);
epicsShareFuncSDDS extern int32_t SDDS_GetRowLimit(void);
epicsShareFuncSDDS extern int32_t SDDS_SetRowLimit(int32_t limit);
epicsShareFuncSDDS extern int32_t SDDS_GotoPage(SDDS_DATASET *SDDS_dataset,int32_t page_number);
epicsShareFuncSDDS extern int32_t SDDS_CheckEndOfFile(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ReadPage(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ReadPageSparse(SDDS_DATASET *SDDS_dataset, uint32_t mode,
                                int32_t sparse_interval,
                                int32_t sparse_offset);
#define SDDS_ReadTable(a) SDDS_ReadPage(a)
epicsShareFuncSDDS extern int32_t SDDS_ReadAsciiPage(SDDS_DATASET *SDDS_dataset, int32_t sparse_interval,
                               int32_t sparse_offset);
epicsShareFuncSDDS extern int32_t SDDS_ReadRecoveryPossible(void);

epicsShareFuncSDDS extern int32_t SDDS_SetColumnFlags(SDDS_DATASET *SDDS_dataset, int32_t column_flag_value);
epicsShareFuncSDDS extern int32_t SDDS_SetRowFlags(SDDS_DATASET *SDDS_dataset, int32_t row_flag_value);
epicsShareFuncSDDS extern int32_t SDDS_GetRowFlag(SDDS_DATASET *SDDS_dataset, int32_t row);
epicsShareFuncSDDS extern int32_t SDDS_GetRowFlags(SDDS_DATASET *SDDS_dataset, int32_t *flag, int32_t rows);
epicsShareFuncSDDS extern int32_t SDDS_BufferedRead(void *target, size_t targetSize, FILE *fp, SDDS_FILEBUFFER *fBuffer);
#if defined(zLib)
epicsShareFuncSDDS extern int32_t SDDS_GZipBufferedRead(void *target, size_t targetSize, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
#endif

/* logic flags for SDDS_AssertRowFlags and SDDS_AssertColumnFlags */
#define SDDS_FLAG_ARRAY  0x001UL
#define SDDS_INDEX_LIMITS 0x002UL
epicsShareFuncSDDS extern int32_t SDDS_AssertRowFlags(SDDS_DATASET *SDDS_dataset, uint32_t mode, ...);
/* modes for SDDS_SetColumnsOfInterest and SDDS_SetRowsOfInterest: */
#define SDDS_NAME_ARRAY 1
#define SDDS_NAMES_STRING 2
#define SDDS_NAME_STRINGS 3
#define SDDS_MATCH_STRING 4
#define SDDS_MATCH_EXCLUDE_STRING 5
#define SDDS_CI_NAME_ARRAY 6
#define SDDS_CI_NAMES_STRING 7
#define SDDS_CI_NAME_STRINGS 8
#define SDDS_CI_MATCH_STRING 9
#define SDDS_CI_MATCH_EXCLUDE_STRING 10

/* logic flags for SDDS_SetColumnsOfInterest, SDDS_SetRowsOfInterest, and SDDS_MatchRowsOfInterest: */
#define SDDS_AND               0x0001UL
#define SDDS_OR                0x0002UL
#define SDDS_NEGATE_MATCH      0x0004UL
#define SDDS_NEGATE_PREVIOUS   0x0008UL
#define SDDS_NEGATE_EXPRESSION 0x0010UL
#define SDDS_INDIRECT_MATCH    0x0020UL
#define SDDS_1_PREVIOUS        0x0040UL
#define SDDS_0_PREVIOUS        0x0080UL
/* used by MatchRowsOfInterest only at this point: */
#define SDDS_NOCASE_COMPARE    0x0100UL

epicsShareFuncSDDS extern int32_t SDDS_MatchColumns(SDDS_DATASET *SDDS_dataset, char ***match, int32_t matchMode, int32_t typeMode, ... );
epicsShareFuncSDDS extern int32_t SDDS_MatchParameters(SDDS_DATASET *SDDS_dataset, char ***match, int32_t matchMode, int32_t typeMode, ... );
epicsShareFuncSDDS extern int32_t SDDS_MatchArrays(SDDS_DATASET *SDDS_dataset, char ***match, int32_t matchMode, int32_t typeMode, ... );
epicsShareFuncSDDS extern int32_t SDDS_Logic(int32_t previous, int32_t match, uint32_t logic);
epicsShareFuncSDDS extern int32_t SDDS_SetColumnsOfInterest(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_AssertColumnFlags(SDDS_DATASET *SDDS_dataset, uint32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetRowsOfInterest(SDDS_DATASET *SDDS_dataset, char *selection_column, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_MatchRowsOfInterest(SDDS_DATASET *SDDS_dataset, char *selection_column, char *label_to_match, int32_t logic);
epicsShareFuncSDDS extern int32_t SDDS_DeleteColumn(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern int32_t SDDS_DeleteParameter(SDDS_DATASET *SDDS_dataset, char *parameter_name);
epicsShareFuncSDDS extern int32_t SDDS_DeleteUnsetColumns(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_CountColumnsOfInterest(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ColumnIsOfInterest(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_ColumnCount(SDDS_DATASET *dataset);
epicsShareFuncSDDS extern int32_t SDDS_ParameterCount(SDDS_DATASET *dataset);
epicsShareFuncSDDS extern int32_t SDDS_ArrayCount(SDDS_DATASET *dataset);
epicsShareFuncSDDS extern int32_t SDDS_CountRowsOfInterest(SDDS_DATASET *SDDS_dataset);
#define SDDS_RowCount(SDDS_dataset) ((SDDS_dataset)->n_rows)
epicsShareFuncSDDS extern int32_t SDDS_DeleteUnsetRows(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_FilterRowsOfInterest(SDDS_DATASET *SDDS_dataset, char *filter_column, double lower, double upper, int32_t logic);
epicsShareFuncSDDS extern int32_t SDDS_ItemInsideWindow(void *data, int32_t index, int32_t type, double lower_limit, double upper_limit);
epicsShareFuncSDDS extern int32_t SDDS_FilterRowsByNumScan(SDDS_DATASET *SDDS_dataset, char *filter_column, uint32_t mode);
#define NUMSCANFILTER_INVERT 0x0001UL
#define NUMSCANFILTER_STRICT 0x0002UL

epicsShareFuncSDDS extern void *SDDS_GetColumn(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern void *SDDS_GetInternalColumn(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern double *SDDS_GetColumnInDoubles(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern float *SDDS_GetColumnInFloats(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern int32_t *SDDS_GetColumnInLong(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern short *SDDS_GetColumnInShort(SDDS_DATASET *SDDS_dataset, char *column_name);
epicsShareFuncSDDS extern void *SDDS_GetNumericColumn(SDDS_DATASET *SDDS_dataset, char *column_name, int32_t desiredType);
epicsShareFuncSDDS extern void *SDDS_GetRow(SDDS_DATASET *SDDS_dataset, int32_t srow_index, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetValue(SDDS_DATASET *SDDS_dataset, char *column_name, int32_t srow_index, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetValueByIndex(SDDS_DATASET *SDDS_dataset, int32_t column_index, int32_t srow_index, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetValueByAbsIndex(SDDS_DATASET *SDDS_dataset, int32_t column_index, int32_t srow_index, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetParameter(SDDS_DATASET *SDDS_dataset, char *parameter_name, void *memory);
epicsShareFuncSDDS extern void *SDDS_GetParameterByIndex(SDDS_DATASET *SDDS_dataset, int32_t index, void *memory);
epicsShareFuncSDDS extern double *SDDS_GetParameterAsDouble(SDDS_DATASET *SDDS_dataset, char *parameter_name, double *data);
epicsShareFuncSDDS extern int32_t *SDDS_GetParameterAsLong(SDDS_DATASET *SDDS_dataset, char *parameter_name, int32_t *data);
epicsShareFuncSDDS extern char *SDDS_GetParameterAsString(SDDS_DATASET *SDDS_dataset, char *parameter_name, char **memory);
epicsShareFuncSDDS extern int32_t SDDS_GetParameters(SDDS_DATASET *SDDS_dataset, ...);
epicsShareFuncSDDS extern void *SDDS_GetFixedValueParameter(SDDS_DATASET *SDDS_dataset, char *parameter_name, void *memory);
epicsShareFuncSDDS extern int32_t SDDS_GetDescription(SDDS_DATASET *SDDS_dataset, char **text, char **contents);

epicsShareFuncSDDS extern void *SDDS_GetMatrixOfRows(SDDS_DATASET *SDDS_dataset, int32_t *n_rows);
epicsShareFuncSDDS extern void *SDDS_GetCastMatrixOfRows(SDDS_DATASET *SDDS_dataset, int32_t *n_rows, int32_t sddsType);
#define SDDS_ROW_MAJOR_DATA 1
#define SDDS_COLUMN_MAJOR_DATA 2
epicsShareFuncSDDS extern void *SDDS_GetMatrixFromColumn(SDDS_DATASET *SDDS_dataset, char *column_name, int32_t dimension1, int32_t dimension2, int32_t mode);
epicsShareFuncSDDS extern void *SDDS_GetDoubleMatrixFromColumn(SDDS_DATASET *SDDS_dataset, char *column_name, int32_t dimension1, int32_t dimension2, int32_t mode);

epicsShareFuncSDDS extern SDDS_ARRAY *SDDS_GetArray(SDDS_DATASET *SDDS_dataset, char *array_name, SDDS_ARRAY *memory);
#define SDDS_POINTER_ARRAY 1
#define SDDS_CONTIGUOUS_DATA 2
epicsShareFuncSDDS extern double *SDDS_GetArrayInDoubles(SDDS_DATASET *SDDS_dataset, char *array_name, int32_t *values);
// epicsShareFuncSDDS extern int32_t SDDS_SetArrayVararg(SDDS_DATASET *SDDS_dataset, char *array_name, int32_t mode, void *data_pointer, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetArrayVararg(SDDS_DATASET *SDDS_dataset, const char *array_name, int32_t mode, void *data_pointer, ...);
epicsShareFuncSDDS extern int32_t SDDS_SetArray(SDDS_DATASET *SDDS_dataset, char *array_name, int32_t mode, void *data_pointer, int32_t *dimension);
epicsShareFuncSDDS extern int32_t SDDS_AppendToArrayVararg(SDDS_DATASET *SDDS_dataset, char *array_name, int32_t mode, void *data_pointer, int32_t elements, ...);


/* error-handling and utility routines: */
epicsShareFuncSDDS extern void *SDDS_Realloc(void *old_ptr, size_t new_size);
epicsShareFuncSDDS extern void *SDDS_Malloc(size_t size);
epicsShareFuncSDDS extern void SDDS_Free(void *mem);
epicsShareFuncSDDS extern void *SDDS_Calloc(size_t nelem, size_t elem_size);
epicsShareFuncSDDS extern int32_t SDDS_NumberOfErrors(void);
epicsShareFuncSDDS extern void SDDS_ClearErrors(void);
epicsShareFuncSDDS extern void SDDS_SetError(char *error_text);
epicsShareFuncSDDS extern void SDDS_Bomb(char *message);
epicsShareFuncSDDS extern void SDDS_Warning(char *message);
epicsShareFuncSDDS extern void SDDS_RegisterProgramName(const char *name);
#define SDDS_VERBOSE_PrintErrors 1
#define SDDS_EXIT_PrintErrors 2
epicsShareFuncSDDS extern void SDDS_PrintErrors(FILE *fp, int32_t mode);
#define SDDS_LAST_GetErrorMessages 0
#define SDDS_ALL_GetErrorMessages 1
epicsShareFuncSDDS extern char **SDDS_GetErrorMessages(int32_t *number, int32_t mode);

epicsShareFuncSDDS extern char **SDDS_GetColumnNames(SDDS_DATASET *SDDS_dataset, int32_t *number);
epicsShareFuncSDDS extern char **SDDS_GetParameterNames(SDDS_DATASET *SDDS_dataset, int32_t *number);
epicsShareFuncSDDS extern char **SDDS_GetAssociateNames(SDDS_DATASET *SDDS_dataset, int32_t *number);
epicsShareFuncSDDS extern char **SDDS_GetArrayNames(SDDS_DATASET *SDDS_dataset, int32_t *number);

epicsShareFuncSDDS extern COLUMN_DEFINITION *SDDS_GetColumnDefinition(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern COLUMN_DEFINITION *SDDS_CopyColumnDefinition(COLUMN_DEFINITION **target, COLUMN_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_FreeColumnDefinition(COLUMN_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_TransferColumnDefinition(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_DefineColumnLikeParameter(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_DefineColumnLikeArray(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_TransferAllColumnDefinitions(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source,
                                              uint32_t mode);


epicsShareFuncSDDS extern PARAMETER_DEFINITION *SDDS_GetParameterDefinition(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern PARAMETER_DEFINITION *SDDS_CopyParameterDefinition(PARAMETER_DEFINITION **target, PARAMETER_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_FreeParameterDefinition(PARAMETER_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_TransferParameterDefinition(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_DefineParameterLikeColumn(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_DefineParameterLikeArray(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
#define SDDS_TRANSFER_KEEPOLD 0x01UL
#define SDDS_TRANSFER_OVERWRITE 0x02UL
epicsShareFuncSDDS extern int32_t SDDS_TransferAllParameterDefinitions(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source, 
                                                 uint32_t mode);

epicsShareFuncSDDS extern ARRAY_DEFINITION *SDDS_GetArrayDefinition(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern ARRAY_DEFINITION *SDDS_CopyArrayDefinition(ARRAY_DEFINITION **target, ARRAY_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_FreeArrayDefinition(ARRAY_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_TransferArrayDefinition(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);
epicsShareFuncSDDS extern int32_t SDDS_TransferAllArrayDefinitions(SDDS_DATASET *SDDS_target, SDDS_DATASET *SDDS_source,
                                             uint32_t mode);


epicsShareFuncSDDS extern ASSOCIATE_DEFINITION *SDDS_GetAssociateDefinition(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern ASSOCIATE_DEFINITION *SDDS_CopyAssociateDefinition(ASSOCIATE_DEFINITION **target, ASSOCIATE_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_FreeAssociateDefinition(ASSOCIATE_DEFINITION *source);
epicsShareFuncSDDS extern int32_t SDDS_TransferAssociateDefinition(SDDS_DATASET *target, SDDS_DATASET *source, char *name, char *newName);

epicsShareFuncSDDS extern int32_t SDDS_GetColumnIndex(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetParameterIndex(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetArrayIndex(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetAssociateIndex(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetColumnType(SDDS_DATASET *SDDS_dataset, int32_t index);
epicsShareFuncSDDS extern int32_t SDDS_GetNamedColumnType(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetParameterType(SDDS_DATASET *SDDS_dataset, int32_t index);
epicsShareFuncSDDS extern int32_t SDDS_GetNamedParameterType(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetArrayType(SDDS_DATASET *SDDS_dataset, int32_t index);
epicsShareFuncSDDS extern int32_t SDDS_GetNamedArrayType(SDDS_DATASET *SDDS_dataset, char *name);
epicsShareFuncSDDS extern int32_t SDDS_GetTypeSize(int32_t type);
epicsShareFuncSDDS extern char *SDDS_GetTypeName(int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_IdentifyType(char *typeName);

#define FIND_ANY_TYPE       0
#define FIND_SPECIFIED_TYPE 1
#define FIND_NUMERIC_TYPE   2
#define FIND_INTEGER_TYPE   3
#define FIND_FLOATING_TYPE  4

epicsShareFuncSDDS extern char *SDDS_FindColumn(SDDS_DATASET *SDDS_dataset,  int32_t mode, ...);
epicsShareFuncSDDS extern char *SDDS_FindParameter(SDDS_DATASET *SDDS_dataset,  int32_t mode, ...);
epicsShareFuncSDDS extern char *SDDS_FindArray(SDDS_DATASET *SDDS_dataset,  int32_t mode, ...);

epicsShareFuncSDDS extern int32_t SDDS_CheckColumn(SDDS_DATASET *SDDS_dataset, char *name, char *units, int32_t type, FILE *fp_message);
epicsShareFuncSDDS extern int32_t SDDS_CheckParameter(SDDS_DATASET *SDDS_dataset, char *name, char *units, int32_t type, FILE *fp_message);
epicsShareFuncSDDS extern int32_t SDDS_CheckArray(SDDS_DATASET *SDDS_dataset, char *name, char *units, int32_t type, FILE *fp_message);
epicsShareFuncSDDS extern int32_t SDDS_VerifyArrayExists(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_VerifyColumnExists(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_VerifyParameterExists(SDDS_DATASET *SDDS_dataset, int32_t mode, ...);
epicsShareFuncSDDS extern int32_t SDDS_PrintCheckText(FILE *fp, char *name, char *units, int32_t type, char *class_name, int32_t error_code);
#define SDDS_CHECK_OKAY 0
#define SDDS_CHECK_OK SDDS_CHECK_OKAY
#define SDDS_CHECK_NONEXISTENT 1
#define SDDS_CHECK_WRONGTYPE  2
#define SDDS_CHECK_WRONGUNITS  3

epicsShareFuncSDDS extern int32_t SDDS_IsActive(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ForceInactive(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_LockFile(FILE *fp, char *filename, char *callerName);
epicsShareFuncSDDS extern int32_t SDDS_FileIsLocked(char *filename);
epicsShareFuncSDDS extern int32_t SDDS_BreakIntoLockedFile(char *filename);

epicsShareFuncSDDS extern int32_t SDDS_CopyString(char **target, char *source);
epicsShareFuncSDDS extern int32_t SDDS_CopyStringArray(char **target, char **source, int32_t n_strings);
epicsShareFuncSDDS extern int32_t SDDS_FreeStringArray(char **string, int32_t strings);
epicsShareFuncSDDS extern int32_t SDDS_VerifyPrintfFormat(char *format_string, int32_t type);
epicsShareFuncSDDS extern int32_t SDDS_HasWhitespace(char *string);
epicsShareFuncSDDS extern char *fgetsSkipComments(char *s, int32_t slen, FILE *fp, char skip_char);
epicsShareFuncSDDS extern char *fgetsSkipCommentsResize(char **s, int32_t *slen, FILE *fp, char skip_char);
#if defined(zLib)
epicsShareFuncSDDS extern char *fgetsGZipSkipComments(char *s, int32_t slen, gzFile *gzfp, char skip_char);
epicsShareFuncSDDS extern char *fgetsGZipSkipCommentsResize(char **s, int32_t *slen, gzFile *gzfp, char skip_char);
#endif
epicsShareFuncSDDS extern void SDDS_CutOutComments(char *s, char cc);
epicsShareFuncSDDS extern void SDDS_EscapeNewlines(char *s);
epicsShareFuncSDDS extern void SDDS_EscapeQuotes(char *s, char quote_char);
epicsShareFuncSDDS extern void SDDS_UnescapeQuotes(char *s, char quote_char);
epicsShareFuncSDDS extern int32_t SDDS_IsQuoted(char *string, char *position, char quotation_mark);
epicsShareFuncSDDS extern int32_t SDDS_GetToken(char *s, char *buffer, int32_t buflen);
epicsShareFuncSDDS extern int32_t SDDS_GetToken2(char *s, char **st, int32_t *strlength, char *buffer, int32_t buflen);
epicsShareFuncSDDS extern int32_t SDDS_PadToLength(char *string, int32_t length);
epicsShareFuncSDDS extern void SDDS_EscapeCommentCharacters(char *string, char cc);
epicsShareFuncSDDS extern void SDDS_InterpretEscapes(char *s);

epicsShareFuncSDDS extern int32_t SDDS_ZeroMemory(void *mem, int32_t n_bytes);
epicsShareFuncSDDS extern int32_t SDDS_SetMemory(void *mem, int32_t n_elements, int32_t data_type, ...);
#define SDDS_PRINT_NOQUOTES 0x0001UL
epicsShareFuncSDDS extern int32_t SDDS_SprintTypedValue(void *data, int32_t index, int32_t type, char *format, char *buffer, uint32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_PrintTypedValue(void *data, int32_t index, int32_t type, char *format, FILE *fp, uint32_t mode);
epicsShareFuncSDDS extern int32_t SDDS_WriteTypedValue(void *data, int32_t index, int32_t type, char *format, FILE *fp);
epicsShareFuncSDDS extern void *SDDS_CastValue(void *data, int32_t index, int32_t data_type, int32_t desired_type, void *memory);
epicsShareFuncSDDS extern void SDDS_RemovePadding(char *s);
epicsShareFuncSDDS extern int32_t SDDS_StringIsBlank(char *s);
epicsShareFuncSDDS extern void *SDDS_AllocateMatrix(int32_t size, int32_t dim1, int32_t dim2);
epicsShareFuncSDDS extern void SDDS_FreeMatrix(void **ptr, int32_t dim1);
epicsShareFuncSDDS extern void SDDS_FreeArray(SDDS_ARRAY *array);
epicsShareFuncSDDS extern void *SDDS_MakePointerArray(void *data, int32_t type, int32_t dimensions, int32_t *dimension);
epicsShareFuncSDDS extern int32_t SDDS_ApplyFactorToParameter(SDDS_DATASET *SDDS_dataset, char *name, double factor);
epicsShareFuncSDDS extern int32_t SDDS_ApplyFactorToColumn(SDDS_DATASET *SDDS_dataset, char *name, double factor);
epicsShareFuncSDDS extern int32_t SDDS_DeleteParameterFixedValues(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_SetDataMode(SDDS_DATASET *SDDS_dataset, int32_t newmode);
epicsShareFuncSDDS extern int32_t SDDS_CheckDataset(SDDS_DATASET *SDDS_dataset, const char *caller);
epicsShareFuncSDDS extern int32_t SDDS_CheckTabularData(SDDS_DATASET *SDDS_dataset, const char *caller);
epicsShareFuncSDDS extern int32_t SDDS_CheckDatasetStructureSize(int32_t size);
#define SDDS_CheckTableStructureSize(a) SDDS_CheckDatasetStructureSize(a)

#define TABULAR_DATA_CHECKS 0x0001UL
epicsShareFuncSDDS uint32_t SDDS_SetAutoCheckMode(uint32_t newMode);

epicsShareFuncSDDS extern int32_t SDDS_FlushBuffer(FILE *fp, SDDS_FILEBUFFER *fBuffer);
epicsShareFuncSDDS extern int32_t SDDS_BufferedWrite(void *target, size_t targetSize, FILE *fp, SDDS_FILEBUFFER *fBuffer);

epicsShareFuncSDDS extern int32_t SDDS_ScanData(char *string, int32_t type, int32_t field_length, void *data, int32_t index, int32_t is_parameter);
epicsShareFuncSDDS extern int32_t SDDS_ScanData2(char *string, char **pstring, int32_t *strlength, int32_t type, int32_t field_length, void *data, int32_t index, int32_t is_parameter);

epicsShareFuncSDDS extern double SDDS_ConvertToDouble(int32_t type, void *data, int32_t index);
epicsShareFuncSDDS extern int32_t SDDS_ConvertToLong(int32_t type, void *data, int32_t index);

epicsShareFuncSDDS extern int32_t SDDS_WriteBinaryString(char *string, FILE *fp, SDDS_FILEBUFFER *fBuffer);
#if defined(zLib)
epicsShareFuncSDDS extern int32_t SDDS_GZipWriteBinaryString(char *string, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
epicsShareFuncSDDS extern int32_t SDDS_GZipBufferedRead(void *target, size_t targetSize, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
epicsShareFuncSDDS extern int32_t SDDS_GZipFlushBuffer(gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
epicsShareFuncSDDS extern int32_t SDDS_GZipBufferedWrite(void *target, size_t targetSize, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
#endif


epicsShareFuncSDDS extern int32_t SDDS_CreateRpnMemory(char *name, short is_string);
epicsShareFuncSDDS extern int32_t SDDS_CreateRpnArray(char *name);

#if defined(RPN_SUPPORT)
epicsShareFuncSDDS extern int32_t SDDS_FilterRowsWithRpnTest(SDDS_DATASET *SDDS_dataset, char *rpn_test);
epicsShareFuncSDDS extern int32_t SDDS_StoreParametersInRpnMemories(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_StoreRowInRpnMemories(SDDS_DATASET *SDDS_dataset, int32_t row);
epicsShareFuncSDDS extern int32_t SDDS_StoreColumnsInRpnArrays(SDDS_DATASET *SDDS_dataset);
epicsShareFuncSDDS extern int32_t SDDS_ComputeColumn(SDDS_DATASET *SDDS_dataset, int32_t column, char *equation);
epicsShareFuncSDDS extern int32_t SDDS_ComputeParameter(SDDS_DATASET *SDDS_dataset, int32_t column, char *equation);
#endif

#define SDDS_BIGENDIAN_SEEN      0x0001UL
#define SDDS_LITTLEENDIAN_SEEN   0x0002UL
#define SDDS_FIXED_ROWCOUNT_SEEN 0x0004UL
#define SDDS_BIGENDIAN         SDDS_BIGENDIAN_SEEN
#define SDDS_LITTLEENDIAN      SDDS_LITTLEENDIAN_SEEN
#define SDDS_FIXED_ROWCOUNT    SDDS_FIXED_ROWCOUNT_SEEN
epicsShareFuncSDDS extern int32_t SDDS_IsBigEndianMachine(void);
void SDDS_SwapShort(short *data);
void SDDS_SwapUShort(unsigned short *data);
epicsShareFuncSDDS extern void SDDS_SwapLong(int32_t *data);
epicsShareFuncSDDS extern void SDDS_SwapULong(uint32_t *data);
void SDDS_SwapFloat(float *data);
void SDDS_SwapDouble(double *data);
epicsShareFuncSDDS extern int32_t SDDS_SwapEndsArrayData(SDDS_DATASET *SDDSin);
epicsShareFuncSDDS extern int32_t SDDS_SwapEndsParameterData(SDDS_DATASET *SDDSin) ;
epicsShareFuncSDDS extern int32_t SDDS_SwapEndsColumnData(SDDS_DATASET *SDDSin);






epicsShareFuncSDDS extern int32_t SDDS_ReadNonNativePage(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_ReadNonNativePageSparse(SDDS_DATASET *SDDS_dataset, uint32_t mode,
                         int32_t sparse_interval,
                         int32_t sparse_offset);
int32_t SDDS_ReadNonNativeBinaryPage(SDDS_DATASET *SDDS_dataset, int32_t sparse_interval, int32_t sparse_offset);
int32_t SDDS_ReadNonNativeBinaryParameters(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_ReadNonNativeBinaryArrays(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_ReadNonNativeBinaryRow(SDDS_DATASET *SDDS_dataset, int32_t row, int32_t skip);
char *SDDS_ReadNonNativeBinaryString(FILE *fp, SDDS_FILEBUFFER *fBuffer, int32_t skip);
#if defined(zLib)
char *SDDS_ReadNonNativeGZipBinaryString(gzFile *gzfp, SDDS_FILEBUFFER *fBuffer, int32_t skip);
#endif



epicsShareFuncSDDS extern int32_t SDDS_WriteNonNativeBinaryPage(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_WriteNonNativeBinaryParameters(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_WriteNonNativeBinaryArrays(SDDS_DATASET *SDDS_dataset);
int32_t SDDS_WriteNonNativeBinaryRow(SDDS_DATASET *SDDS_dataset, int32_t row);

int32_t SDDS_WriteNonNativeBinaryString(char *string, FILE *fp, SDDS_FILEBUFFER *fBuffer);
#if defined(zLib)
int32_t SDDS_GZipWriteNonNativeBinaryString(char *string, gzFile *gzfp, SDDS_FILEBUFFER *fBuffer);
#endif

#if defined(SDDS_MPI_IO) && SDDS_MPI_IO
  /* SDDSmpi_output.c */
  char *BlankToNull(char *string);
  epicsShareFuncSDDS extern void SDDS_MPI_BOMB(char *text, MPI_File *mpi_file);
  void SDDS_MPI_GOTO_ERROR(FILE *fp, char *str, int32_t mpierror, int32_t exit);
  epicsShareFuncSDDS extern int32_t SDDS_MPI_File_Open(MPI_DATASET *MPI_dataset, char *filename, unsigned long flags);
  char *SDDS_CreateNamelistField(char *name, char *value);
  char *SDDS_CreateDescription(char *text, char *contents);
  char *SDDS_CreateParameterDefinition(PARAMETER_DEFINITION *parameter_definition);
  char *SDDS_CreateColumnDefinition(COLUMN_DEFINITION *column_definition);
  char *SDDS_CreateArrayDefinition(ARRAY_DEFINITION *array_definition);
  char *SDDS_CreateAssociateDefinition(ASSOCIATE_DEFINITION *associate_definition);
  char *SDDS_CreateDataMode(DATA_MODE *data_mode);
#define SDDS_MPI_WriteTable(a) SDDS_MPI_WritePage(a)
  epicsShareFuncSDDS extern int32_t SDDS_MPI_WriteLayout(SDDS_DATASET *MPI_dataset);
  epicsShareFuncSDDS extern int32_t SDDS_MPI_WritePage(SDDS_DATASET *MPI_dataset);
  MPI_Datatype Convert_SDDStype_To_MPItype(int32_t SDDS_type);
  epicsShareFuncSDDS extern int32_t SDDS_MPI_Terminate(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_InitializeOutput(SDDS_DATASET *MPI_dataset, char *description, char *contents, char *filename, unsigned long flags, short column_major);
  int32_t SDDS_MPI_InitializeCopy(SDDS_DATASET *MPI_dataset_target, SDDS_DATASET *SDDS_source, char *filename, short column_major);
  
  /*SDDS_MPI_binary.c writing data*/
  int32_t SDDS_CheckStringTruncated(void);
  void SDDS_StringTuncated(void);
  int32_t SDDS_SetDefaultStringLength(int32_t newValue);
  int32_t SDDS_MPI_WriteBinaryPage(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_WriteBinaryString(SDDS_DATASET *MPI_dataset, char *string);
  int32_t SDDS_MPI_WriteBinaryParameters(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_WriteBinaryArrays(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_WriteBinaryRow(SDDS_DATASET *MPI_dataset, int32_t row);
   int32_t SDDS_MPI_WriteNonNativeBinaryPage(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_WriteNonNativeBinaryString(SDDS_DATASET *MPI_dataset, char *string);
  int32_t SDDS_MPI_WriteNonNativeBinaryParameters(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_WriteNonNativeBinaryArrays(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_WriteNonNativeBinaryRow(SDDS_DATASET *MPI_dataset, int32_t row);
  MPI_Offset SDDS_MPI_Get_Column_Size(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_CollectiveWriteByRow(SDDS_DATASET *SDDS_dataset);
  int32_t SDDS_MPI_Get_Title_Size(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_BufferedWrite(void *target, size_t targetSize, SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_FlushBuffer(SDDS_DATASET *MPI_Dataset);
  int32_t SDDS_MPI_GetTotalRows(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_CountRowsOfInterest(SDDS_DATASET *SDDS_dataset, int32_t start_row, int32_t end_row);
  int32_t SDDS_MPI_WriteContinuousBinaryPage(SDDS_DATASET *MPI_dataset);
  MPI_Offset SDDS_MPI_GetTitleOffset(SDDS_DATASET *MPI_dataset);
  /*SDDS_MPI_binary.c reading data*/
  int32_t SDDS_MPI_BufferedRead(void *target, size_t targetSize, SDDS_DATASET *MPI_dataset, SDDS_FILEBUFFER *fBuffer);
  int32_t SDDS_MPI_ReadBinaryPage(SDDS_DATASET *MPI_dataset);
  char *SDDS_MPI_ReadNonNativeBinaryString(SDDS_DATASET *MPI_dataset, SDDS_FILEBUFFER *fBuffer, int32_t skip);
  int32_t SDDS_MPI_ReadBinaryParameters(SDDS_DATASET *MPI_dataset, SDDS_FILEBUFFER *fBuffer);
  int32_t SDDS_MPI_ReadBinaryArrays(SDDS_DATASET *MPI_dataset, SDDS_FILEBUFFER *fBuffer);
  int32_t SDDS_MPI_ReadBinaryRow(SDDS_DATASET *MPI_dataset, int32_t row, int32_t skip);
  int32_t SDDS_MPI_ReadNonNativeBinaryParameters(SDDS_DATASET *SDDS_dataset);
  int32_t SDDS_MPI_ReadNonNativeBinaryArrays(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_ReadNonNativeBinaryRow(SDDS_DATASET *MPI_dataset, int32_t row, int32_t skip);
  int32_t SDDS_MPI_ReadBinaryPage(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_ReadNonNativePage(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_ReadNonNativePageSparse(SDDS_DATASET *MPI_dataset, uint32_t mode);
  int32_t SDDS_MPI_ReadNonNativeBinaryPage(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_MPI_ReadNonNativeBinaryRow(SDDS_DATASET *MPI_dataset, int32_t row, int32_t skip);
  int32_t SDDS_MPI_BufferedReadBinaryTitle(SDDS_DATASET *MPI_dataset);
  int32_t SDDS_SetDefaultTitleBufferSize(int32_t newSize);
  int32_t SDDS_MPI_WriteBinaryPageByColumn(SDDS_DATASET *MPI_dataset);
  epicsShareFuncSDDS extern void SDDS_MPI_Setup(SDDS_DATASET *SDDS_dataset, int32_t parallel_io, int32_t n_processors, int32_t myid, MPI_Comm comm, short master_read);
  
  /*SDDSmpi_input.c */
  epicsShareFuncSDDS extern int32_t SDDS_MPI_ReadPage(SDDS_DATASET *MPI_dataset);
  epicsShareFuncSDDS extern int32_t SDDS_MPI_InitializeInput(SDDS_DATASET *MPI_dataset, char *filename);
  epicsShareFuncSDDS extern int32_t SDDS_MPI_InitializeInputFromSearchPath(SDDS_DATASET *MPI_dataset, char *file);
 #define SDDS_MPI_TotalRowCount(SDDS_DATASET) ((SDDS_DATASET)->MPI_dataset->total_rows) 
#endif

#ifdef __cplusplus
}
#endif

#endif
