#include <time.h>
#include "bstat.h"
#include "statistic.h"

byte    (*_entropy)(void) = 0ul;

void entropy_f(byte *src, qword len) {
    assert(_entropy != NULL);

    for (qword i = 0; i < len; ++i)
        src[i] = _entropy();
}

qword entropy_s(byte *src, qword len){
    assert(_entropy != NULL);

    qword sum = 0;
    for (qword i = 0; i < len; ++i) {
        src[i] = _entropy();
        sum += src[i];
    }
    return sum;
}

void entropy_m(byte *src, qword len, byte mod){
    assert(_entropy != NULL);

    for (qword i = 0; i < len; ++i)
        src[i] = _entropy() % mod;
}



byte ch(double *arr, byte len){
    assert(_entropy != NULL);

    double val = (double)_entropy() / UINT8_MAX;
    double sum = 0;
    for (byte i = 0; i < len; ++i) {
        sum += arr[i];
        if(val <= sum)
            return i;
    }
    return (len - 1);
}

void checkRAM(sequence *ctx){
    if(ctx == NULL)
        assert(!"memory error");
    struct sysinfo info;
    int rv =sysinfo(&info);
    if(ctx == NULL || (ctx->T * 8) > info.freeram)
        assert(!"memory error");
}

void init(){
    srand(time(NULL));
    _entropy = rand;
}