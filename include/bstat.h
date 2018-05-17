#ifndef __BSTAT_BSTAT_H
#define __BSTAT_BSTAT_H

#include <stdint.h>
#include <assert.h>
#include <stdio.h>

#ifdef __linux__
    #include <sys/sysinfo.h>
#endif //__linux__

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
    byte  *array;       /**< The pointer to sequence. */
} sequence;

void free_sequence(sequence *ctx);
void init_sequence(sequence *ctx, qword t);
void checkRAM(sequence *ctx);

typedef struct istream_s{
    FILE *file;
    byte m;
    word buffer;
    byte buffer_size;
    byte bits;
} istream;

byte    iopen(istream *ctx, const char* filename, byte m);
void    iclose(istream *ctx);
byte    iopened(istream *ctx);
byte    iend(istream *ctx);
byte    iget(istream *ctx);
qword   isize(istream *ctx);
/////

typedef struct ostream_s{
    FILE *file;
    byte m;
    byte isASCII;
    word buffer;
    byte buffer_size;
    byte bits;
} ostream;

byte    oopen(ostream *ctx, const char* filename, byte isASCII);
byte    osetm(ostream *ctx, byte m);
void    oclose(ostream *ctx);
byte    oopened(ostream *ctx);
void    oput(ostream *ctx, byte in);


byte read_file_2_sequence(sequence *seq, const char *filename, byte m);


extern byte    (*_entropy)(void);


void    init();
void    entropy_f(byte *, qword);
qword   entropy_s(byte *, qword);
void    entropy_m(byte *, qword, byte);

byte    ch(double *, byte);





#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_BSTAT_H
