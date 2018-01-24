#ifndef __BSTAT_BSTAT_H
#define __BSTAT_BSTAT_H

#include <stdint.h>
#include <assert.h>
#include <stdio.h>
//#include <math.h>

#include <stdlib.h>
#ifdef __cplusplus
extern "C"{
#endif



#define SUCCESS         0
#define ERROR           1
#define CRASH_ERROR     2

#define _memcheck(src, size) do { if (((src) = malloc(size)) == NULL) assert(!CRASH_ERROR); } while(0)



typedef unsigned char byte;
typedef uint16_t      word;
typedef uint32_t      dword;
typedef uint64_t      qword;





extern byte    (*_entropy)(void);

void    entropy_f(byte *, qword);
qword   entropy_s(byte *, qword);
void    entropy_m(byte *, qword, byte);

byte    ch(double *, byte);





#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_BSTAT_H
