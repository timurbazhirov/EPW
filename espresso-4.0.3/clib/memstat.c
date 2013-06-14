/*
  Copyright (C) 2002 FPMD group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

#include "c_defs.h"

/* 
  This function return the numer of kilobytes allocated
  by the calling process. 
  Auhor: Carlo Cavazzoni.
*/

void F77_FUNC(memstat,MEMSTAT)(int *kilobytes)
{
#if defined(HAVE_MALLINFO) && !defined(__QK_USER__) && !defined(__SOLARIS) 
#include <malloc.h>
  struct mallinfo info;  
  info = mallinfo();

#if defined(__AIX)
  *kilobytes = (info.arena) / 1024 ;
#else
  *kilobytes = (info.arena + info.hblkhd) / 1024 ;
#endif

#else
  *kilobytes = -1;
#endif
}
