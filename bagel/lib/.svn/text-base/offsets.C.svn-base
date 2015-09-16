/*
 *
 *  Copyright Peter Boyle and Glasgow University 2000.
 *
 *  It is provided under the GNU pubic License V2
 *  It is provided as is and is not guaranteed fit for any purpose.
 *
 */


#include "processor.h"
extern struct processor *PROC;

int SizeofDatum(enum Datum type)
{
  int type_size;

  test_proc();

  switch ( type ) {
  case Integer:
    type_size = PROC->I_size; break;

  case ShortInteger:
    type_size = PROC->IS_size; break;

  case Single:
    type_size = PROC->FSP_size; break;

  case Double:
    type_size = PROC->FP_size; break;

  case Byte:
    type_size = 1; break;

  case Half:
    type_size = 2; break;

  default:
    exit(-1);break;
  }
  return type_size;
}

/* {...IMMEDIATE offset management routines*/

struct offset *offset_list = NULL;
int def_offset ( int units,enum Datum type,const char *str )
{
  struct offset *newo = new offset;
  struct offset *search;
  int handle;

  int type_size = SizeofDatum(type);

  newo -> offst = units * type_size;
  newo -> name = str;
  newo -> next = NULL;

  handle = 0;
  if ( offset_list == NULL ) {
    offset_list = newo;
    return( handle);
  }
  search  = offset_list;
  while ( search-> next != NULL ) {
    handle ++;
    search = search -> next;
  }

  /*Add to end of list*/
  search->next = newo;
  handle ++;

  return(handle);
}
int get_offset( int handle )
{
  if ( handle < 0 ) {
    fprintf(stderr,"Bad handle %d passed to get_offset\n",handle);
    exit(-1);
  }
  return( get_offset_struct(handle) -> offst );
}
int get_offset_handle (int off,enum Datum size)
{
  int handle = 0;

  struct offset *search;

  search  = offset_list;
  while ( search != NULL ) {

    if ( search-> offst == off ) return handle;  // FIXME: assuming size==Byte here?!

    search = search -> next;
    handle ++;

  }

  char * buf = (char *)malloc(80);
  sprintf(buf,"auto_%d",off);

  return def_offset (off,size,buf);
}


const char * get_offset_name ( int handle )
{
  return( get_offset_struct(handle) -> name );
}

struct offset * get_offset_struct( int handle )
{
  struct offset * search;
  int count;

  count = 0 ;
  search = offset_list;

  while ( count != handle && search != NULL) {
    search = search -> next;
    count ++;
  }
  return ( search);
}
/* ...}*/
