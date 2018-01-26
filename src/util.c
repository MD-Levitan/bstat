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

    double val = (double)_entropy() / INT8_MAX;
    double sum = 0;
    for (byte i = 0; i < len; ++i) {
        sum += arr[i];
        if(val <= sum)
            return i;
    }
    return (len - 1);
}

double standart_deviation_matrix(const double **s1, const double **s2, qword dim1, qword dim2){
    double std_dev = 0;
    for (qword i = 0; i < dim1; ++i)
        for (qword j = 0; j < dim2; ++j)
            std_dev += (s1[i][j] - s2[i][j]) * (s1[i][j] - s2[i][j]);
    return std_dev;
}