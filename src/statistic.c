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



byte statistic_chi_square_cmm(byte N, double *n, cmm_model *estimation_cmm, cmm_model *prediction_cmm,
                              double threshold, double *stat){
    if(estimation_cmm == NULL || n == NULL || stat == NULL){
        return ERROR;
    }

    byte need_free = 0;
    if(prediction_cmm == NULL){
        need_free = 1;
        cmm_model model;
        init_cmm_model(&model, estimation_cmm->N);
        generate_uniform_cmm_model(&model);
        prediction_cmm = &model;
    }
    double _stat = 0;
    for(byte i = 0; i < estimation_cmm->N; ++i) {
        for (byte j = 0; j < estimation_cmm->N; ++j) {
            double var = (estimation_cmm->P[i][j] - prediction_cmm->P[i][j]);
            if (prediction_cmm->P[i][j])
                _stat += var * var * n[i] / prediction_cmm->P[i][j];
        }
    }
    *stat = _stat;

    if(need_free)
        prediction_cmm = NULL;
}


byte statistic_chi_square_cmms(byte N, byte S, double *n, cmms_model *estimation_cmms, cmms_model *prediction_cmms,
                               double threshold, double *stat){
    if(estimation_cmms == NULL || n == NULL || stat == NULL){
        return ERROR;
    }

    byte need_free = 0;
    if(prediction_cmms == NULL){
        need_free = 1;
        cmms_model model;
        init_cmms_model(&model, estimation_cmms->N, estimation_cmms->S);
        generate_uniform_cmms_model(&model);
        prediction_cmms = &model;
    }
    double _stat = 0;
    for(byte i = 0; i < pow(estimation_cmms->N, estimation_cmms->S); ++i) {
        for (byte j = 0; j < estimation_cmms->N; ++j) {
            double var = (estimation_cmms->P[i][j] - prediction_cmms->P[i][j]);
            if (prediction_cmms->P[i][j])
                _stat += var * var * n[i] / prediction_cmms->P[i][j];
        }
    }
    *stat = _stat;

    if(need_free)
        prediction_cmms = NULL;
}

/*
 * Statistical hypothesis testing. Use likelihood-ratio statistic.
 * :param sequence: sequence of elements in range (0, A).
 * :param estimation_hmm: params of hidden markov model(estimation).
 * :param prediction_hmm: params of hidden markov model(predictable).
 * :param threshold: threshold probability — the significance level.
 */
byte statistic_likelihood_hmm(sequence *seq,  hmm_model *estimation_hmm, hmm_model *prediction_hmm,
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

//byte statistic_likelihood_dcmm(sequence *seq, dcmm_model *estimation_hmm, dcmm_model *prediction_hmm,
//                               double threshold, double *stat){
//    if(seq == NULL || estimation_hmm == NULL || stat == NULL){
//        return ERROR;
//    }
//
//    byte need_free = 0;
//    if(prediction_hmm == NULL){
//        need_free = 1;
//        dcmm_model model;
//        init_dcmm_model(&model, estimation_hmm->N, estimation_hmm->M);
//        generate_uniform_dcmm_model(&model);
//        prediction_hmm = &model;
//    }
//    double ** alphaset;
//    double * alphaset_v;
//    double est_likehood, pred_likehood;
//
//    init_set(&alphaset, seq, estimation_hmm);
//    init_set_v(&alphaset_v, seq);
//
//    forward_algorithm(seq, estimation_hmm, alphaset, alphaset_v);
//    est_likehood = estimation_sequence_forward(seq, estimation_hmm, alphaset, alphaset_v);
//
//    forward_algorithm(seq, prediction_hmm, alphaset, alphaset_v);
//    pred_likehood = estimation_sequence_forward(seq, prediction_hmm, alphaset, alphaset_v);
//
//    free_set(alphaset, seq);
//    free_set_v(alphaset_v);
//
////    est_likehood = !est_likehood ? 0 : log2(est_likehood);
////    pred_likehood = !pred_likehood ? 0 : log2(pred_likehood);
//    *stat = -2 * (pred_likehood - est_likehood);
//
//    if(need_free)
//        prediction_hmm = NULL;
//
////from scipy.stats import chi2
////        level = chi2.isf(threshold, sequence.N * (sequence.N - 1) + sequence.N * (sequence.M - 1))
////
//
//}