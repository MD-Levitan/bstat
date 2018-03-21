#ifndef __BSTAT_STATISTIC_H
#define __BSTAT_STATISTIC_H

#include "bstat.h"
#include "hmm.h"
#include "cmm.h"

#ifdef __cplusplus
extern "C"{
#endif

double standart_deviation_matrix(const double **s1, const double **s2, qword dim1, qword dim2);
byte statistic_chi_square(int N, double *n, double** estimation_P, double** prediction_P,
                          double threshold, double *stat);
byte statistic_likelihood(struct sequence_hmm_s *seq, struct hmm_s *estimation_hmm, struct hmm_s *prediction_hmm,
                          double threshold, double *stat);




#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_STATISTIC_H
