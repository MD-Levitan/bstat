#ifndef __BSTAT_DCMM_H
#define __BSTAT_DCMM_H

#include "bstat.h"
#include "hmm.h"

#ifdef __cplusplus
extern "C"{
#endif


/**
 * @brief Params of DCMM(Double Chain Markov Model).
 *
 *
 */
typedef struct dcmm_s{
    byte N;             /**< Size of the set  of hidden states. */
    byte M;             /**< Size of the set  of visible states. */
    double *Pi;         /**< Initial vector of hidden states. */
    double **P;         /**< One-step transition matrix of hidden states. */
    double ***C;         /**<A set of transition matrix from hidden and visible states to new visible state. */
} dcmm_model;


void init_dcmm_model(dcmm_model *ctx, byte n, byte m);
void copy_dcmm_model(dcmm_model *dist, dcmm_model *src);
void set_dcmm_model(dcmm_model *dist, double *Pi, double** P, double*** C);
void free_dcmm_model(dcmm_model *model);

/**
 * @brief Generate params of DCMM.
 *
 * @param model The pointer of model.
 * @param param If type != 2 param must be NULL. @p param = {num1, num2, ind}.
 *              Divide domain of a mode into #num1(P) * #num2(C) parts. Generate DCMM in #ind part of domain.
 * @return      Return #Success or #Error.
 *
 */
byte generate_dcmm_model(dcmm_model *model, byte type, uint32_t* param);
#define generate_random_dcmm_model(model) generate_dcmm_model(model, 0, NULL)
#define generate_uniform_dcmm_model(model) generate_dcmm_model(model, 1, NULL)


/**
 *
 * @param seq
 * @param t
 * @param n
 * @param m
 * @param type
 * @return
 */
byte generate_dcmm_sequence(sequence *seq, dcmm_model *model);



byte forward_algorithm_d(sequence *seq, dcmm_model *model, double **set, double *set_v);
byte backward_algorithm_d(sequence *seq, dcmm_model *model, double **set, double *set_v);
double estimation_sequence_forward_d(sequence *seq, dcmm_model *model, double **set, double *set_v);
void double_probability_norm_d(sequence *seq, dcmm_model *model, double estimation_seq, double **alphaset,
                             double *alhaset_v, double **betaset, double *betaset_v, double  ***ksiset);
void marginal_probability_norm_d(sequence *seq, dcmm_model *model, double estimation_seq, double **alphaset,
                                double *alhaset_v, double **betaset, double *betaset_v, double  **gammaset);

void estimation_model_d(sequence *seq, dcmm_model *model, double eps, double *likehood);
void estimation_model_gl_d(sequence *seq, dcmm_model *model, uint32_t iter, double eps, double *likehood); //search for global maximum of likehood.

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_DCMM_H
