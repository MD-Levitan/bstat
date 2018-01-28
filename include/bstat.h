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

#define _check(_) do{ if (_ != SUCCESS) return ERROR; } while(0)
#define _memcheck(src, size) do { if (((src) = malloc(size)) == NULL) assert(!CRASH_ERROR); } while(0)
#define __memcheck(src, num, size) do { if (((src) = calloc(num, size)) == NULL) assert(!CRASH_ERROR); } while(0)


typedef unsigned char byte;
typedef uint16_t      word;
typedef uint32_t      dword;
typedef uint64_t      qword;


typedef struct sequence_s{
    qword T;            /**< The length of sequence. */
    byte  m;            /**< The size of alphabet of seq, if model == NULL. */
    byte* array;        /**< The pointer to sequence. */
} sequence;

void free_sequence(sequence *ctx);
void init_sequence(sequence *ctx, qword t);


extern byte    (*_entropy)(void);

void    entropy_f(byte *, qword);
qword   entropy_s(byte *, qword);
void    entropy_m(byte *, qword, byte);

byte    ch(double *, byte);





#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_BSTAT_H
