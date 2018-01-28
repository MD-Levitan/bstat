#include "statistic.h"

double standart_deviation_matrix(const double **s1, const double **s2, qword dim1, qword dim2){
    double std_dev = 0;
    for (qword i = 0; i < dim1; ++i)
        for (qword j = 0; j < dim2; ++j)
            std_dev += (s1[i][j] - s2[i][j]) * (s1[i][j] - s2[i][j]);
    return std_dev;
}

/*
 * """
Statistical hypothesis testing. Use chi-square statistic.
:param N: size of alphabet.
:param n: vector, with numbers of times that state i observed in sample( i in range (0, N)).
:param estimation_P: one-step transition matrix(estimation).
:param prediction_P: one_step transition matrix(predictable).
:param threshold: threshold probability — the significance level.
"""
 */
//
byte statistic_chi_square(int N, double *n, double** estimation_P, double** prediction_P,
                            double threshold, double *stat);
//
//def statistic_chi_square(N, n, estimation_P, prediction_P=None, threshold=0.05):
//if prediction_P is None:
//        import numpy as np
//        prediction_P = np.full((N, N), 1 / N)
//if len(n) != N:
//return
//
//degree_m = sum((prediction_P != 0).sum(1))
//
//stat = 0.0
//for i in range(0, N):
//stat += sum([(estimation_P[j][i] - prediction_P[j][i]) * (estimation_P[j][i] - prediction_P[j][i]) * n[i]
/// prediction_P[j][i] for j in range(0, N) if prediction_P[j][i] != 0])
//
//from scipy.stats import chi2
//        level = chi2.isf(threshold, degree_m - N)
//
//print("${\gamma}_T = $" + str(stat) + "$, \delta " + str(level))
//
//return [True if stat <= level else False, stat, level]
//


/*
 * Statistical hypothesis testing. Use likelihood-ratio statistic.
 * :param sequence: sequence of elements in range (0, A).
 * :param estimation_hmm: params of hidden markov model(estimation).
 * :param prediction_hmm: params of hidden markov model(predictable).
 * :param threshold: threshold probability — the significance level.
 */
byte statistic_likelihood(hmm_seq *seq, hmm_model *estimation_hmm, hmm_model *prediction_hmm,
                          double threshold, double *stat){
    if(seq == NULL || estimation_hmm == NULL || stat == NULL){
        return ERROR;
    }
    
    byte need_free = 0;
    if(prediction_hmm == NULL){
        need_free = 1;
        hmm_model model;
        init_hmm_model(&model, estimation_hmm->N, estimation_hmm->M);
        generate_uniform_hmm_model(&model);
        prediction_hmm = &model;
    }
    double ** alphaset;
    double * alphaset_v;
    double est_likehood, pred_likehood;
    
    init_set(&alphaset, seq, estimation_hmm);
    init_set_v(&alphaset_v, seq);
    
    forward_algorithm(seq, estimation_hmm, alphaset, alphaset_v);
    est_likehood = estimation_sequence_forward(seq, estimation_hmm, alphaset, alphaset_v);
    
    forward_algorithm(seq, prediction_hmm, alphaset, alphaset_v);
    pred_likehood = estimation_sequence_forward(seq, prediction_hmm, alphaset, alphaset_v);
    
    free_set(alphaset, seq);
    free_set_v(alphaset_v);
    
//    est_likehood = !est_likehood ? 0 : log2(est_likehood);
//    pred_likehood = !pred_likehood ? 0 : log2(pred_likehood);
    *stat = -2 * (pred_likehood - est_likehood);
    
    if(need_free)
        prediction_hmm = NULL;

//from scipy.stats import chi2
//        level = chi2.isf(threshold, sequence.N * (sequence.N - 1) + sequence.N * (sequence.M - 1))
//

}
