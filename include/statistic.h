#ifndef __BSTAT_STATISTIC_H
#define __BSTAT_STATISTIC_H


#include "hmm.h"
#include "cmm.h"
#include "cmms.h"
#include "dcmm.h"

typedef struct hmm_s hmm_model;
typedef struct dcmm_s dcmm_model;

#ifdef __cplusplus
extern "C"{
#endif

double standart_deviation_matrix(const double **s1, const double **s2, qword dim1, qword dim2);

byte statistic_chi_square_cmm(byte N, double *n, cmm_model *estimation_cmm, cmm_model *prediction_cmm,
                          double threshold, double *stat);

byte statistic_chi_square_cmms(byte N, byte S, double *n, cmms_model *estimation_cmms, cmms_model *prediction_cmms,
                               double threshold, double *stat);


byte statistic_likelihood_hmm(sequence *seq,  hmm_model *estimation_hmm, hmm_model *prediction_hmm,
                          double threshold, double *stat);

byte statistic_likelihood_dcmm(sequence *seq, dcmm_model *estimation_hmm, dcmm_model *prediction_hmm,
                          double threshold, double *stat);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_STATISTIC_H
