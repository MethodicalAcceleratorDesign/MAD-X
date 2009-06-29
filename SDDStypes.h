/*************************************************************************\
* Copyright (c) 2002 The University of Chicago, as Operator of Argonne
* National Laboratory.
* Copyright (c) 2002 The Regents of the University of California, as
* Operator of Los Alamos National Laboratory.
* This file is distributed subject to a Software License Agreement found
* in the file LICENSE that is included with this distribution. 
\*************************************************************************/

/* file: SDDS.h
 * purpose: SDDS data types
 *          format
 *
 * Michael Borland, 1993
 $Log: not supported by cvs2svn $
 Revision 1.7  2006/08/31 15:06:55  soliday
 Updated to work with SDDS2

 Revision 1.6  2002/08/14 15:40:13  soliday
 Added Open License

 Revision 1.5  1999/02/19 22:52:26  borland
 Added SDDS_ANY_INTEGER_TYPE for use with SDDS_CheckXXX procedures.

 Revision 1.4  1997/12/19 16:55:46  borland
 Fixed SDDS_RowCount macro (more parentheses).  Added prototype for
 SDDS_Malloc.  Added new "type": SDDS_ANY_FLOATING_TYPE.

 * Revision 1.3  1995/09/06  14:12:01  saunders
 * First test release of SDDS1.5
 *
 */

#if !defined(_SDDSTYPES_)

#define _SDDSTYPES_ 1

#define SDDS_DOUBLE    1
#define SDDS_FLOAT     2
#define SDDS_LONG      3
#define SDDS_ULONG     4
#define SDDS_SHORT     5
#define SDDS_USHORT    6
#define SDDS_STRING    7
#define SDDS_CHARACTER 8
#define SDDS_NUM_TYPES 8
#define SDDS_INTEGER_TYPE(type) ((type)==SDDS_LONG || (type)==SDDS_ULONG || (type)==SDDS_SHORT || (type)==SDDS_USHORT)
#define SDDS_FLOATING_TYPE(type) ((type)==SDDS_DOUBLE || (type)==SDDS_FLOAT)
#define SDDS_NUMERIC_TYPE(type) (SDDS_INTEGER_TYPE(type) || SDDS_FLOATING_TYPE(type))
#define SDDS_VALID_TYPE(type) (type>=1 && type<=SDDS_NUM_TYPES)

/* used by SDDS_Check*() routines  */
#define SDDS_ANY_NUMERIC_TYPE (SDDS_NUM_TYPES+1)
#define SDDS_ANY_FLOATING_TYPE   (SDDS_NUM_TYPES+2)
#define SDDS_ANY_INTEGER_TYPE   (SDDS_NUM_TYPES+3)

#endif
