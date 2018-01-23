
#include "bstat.h"

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

    double val = (double)_entropy() / INT8_MAX;
    double sum = 0;
    for (byte i = 0; i < len; ++i) {
        sum += arr[i];
        if(val <= sum)
            return i;
    }
    return (len - 1);
}
