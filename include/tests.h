#ifndef __BSTAT_TESTS_H
#define __BSTAT_TESTS_H
#include "bstat.h"
#include "cmm.h"
#include "cmms.h"
#include "hmm.h"
#include "dcmm.h"
#include "statistic.h"
#include "time.h"



double *statistic_cmms(byte N, byte S, qword len, byte alg);
double *statistic_cmm(byte N, qword len, byte alg);


double *statistic_cmms_2(byte N, byte S, qword len, byte alg);

double statistic_cmms_3(byte N, byte S, qword len, byte alg);

double statistic_dcmm(byte N, byte M, qword len);
double statistic_dcmm_3(byte N, byte M, qword len);


#endif //__BSTAT_TESTS_H
