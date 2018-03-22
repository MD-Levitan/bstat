#ifndef __BSTAT_HMM_H
#define __BSTAT_HMM_H

#include "bstat.h"
#include "statistic.h"

#ifdef __cplusplus
extern "C"{
#endif


/**
 * @brief Params of HMM(Hidden Markov Model).
 *
 *
 */
typedef struct hmm_s{
    byte N;             /**< Size of the set  of hidden states. */
    byte M;             /**< Size of the set  of visible states. */
    double *Pi;         /**< Initial vector of hidden states. */
    double **P;         /**< One-step transition matrix of hidden states. */
    double **C;         /**< Transition matrix from hidden to visible states. */
} hmm_model;


void init_hmm_model(hmm_model *ctx, byte n, byte m);
void copy_hmm_model(hmm_model *dist, hmm_model *src);
void free_hmm_model(hmm_model *model);

/**
 * @brief Generate params of HMM.
 *
 * @param model The pointer of model.
 * @param param If type != 2 param must be NULL. @p param = {num1, num2, ind}.
 *              Divide domain of a mode into #num1(P) * #num2(C) parts. Generate HMM in #ind part of domain.
 * @return      Return #Success or #Error.
 *
 */
byte generate_hmm_model(hmm_model *model, byte type, uint32_t* param);
#define generate_random_hmm_model(model) generate_hmm_model(model, 0, NULL)
#define generate_uniform_hmm_model(model) generate_hmm_model(model, 1, NULL)

/**
 *
 * @param seq
 * @param t
 * @param n
 * @param m
 * @param type
 * @return
 */
byte generate_hmm_sequence(sequence *seq, hmm_model *model);



byte forward_algorithm(sequence *seq, hmm_model *model, double **set, double *set_v);
byte backward_algorithm(sequence *seq, hmm_model *model, double **set, double *set_v);
double estimation_sequence_forward(sequence *seq, hmm_model *model, double **set, double *set_v);
void double_probability_norm(sequence *seq, hmm_model *model, double estimation_seq, double **alphaset,
                             double *alhaset_v, double **betaset, double *betaset_v, double  ***ksiset);
void marginaol_probability_norm(sequence *seq, hmm_model *model, double estimation_seq, double **alphaset,
                                double *alhaset_v, double **betaset, double *betaset_v, double  **gammaset);

void estimation_model(sequence *seq, hmm_model *model, double eps, double *likehood);
void estimation_model_gl(sequence *seq, hmm_model *model, uint32_t iter, double eps, double *likehood); //search for global maximum of likehood.

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_HMM_H
