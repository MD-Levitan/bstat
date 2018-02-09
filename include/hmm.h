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

typedef struct sequence_hmm_s{
    hmm_model *model;   /**< HMM */
    qword T;            /**< The length of sequence. */
    byte  m;            /**< The size of alphabet of seq, if model == NULL. */
    byte* array;        /**< The pointer to sequence. */
} hmm_seq;

void init_hmm_model(hmm_model *ctx, byte n, byte m);
void init_hmm_seq(hmm_seq *ctx, qword t, hmm_model *ctx2);
void copy_hmm_model(hmm_model *dist, hmm_model *src);
void free_hmm_model(hmm_model *model);
void free_hmm_seq(hmm_seq *seq);

/**
 * @brief Generate params of HMM.
 *
 * @param model The pointer of model.
 * @param type  The type of generation params of HMM.
 *              @p type = {0, 1, 2}, where 0 - random, 1 - uniform, 2 - in section.
 * @param param If type != 2 param must be NULL. @p param = {num1, num2, ind}.
 *              Divide domain of a mode into #num1(P) * #num2(C) parts. Generate HMM in #ind part of domain.
 * @return      Return #Success or #Error.
 *
 */
byte generate_hmm_model(hmm_model *model, byte type, uint32_t* param);
#define generate_random_hmm_model(model) generate_hmm_model(model, 0)
#define generate_uniform_hmm_model(model) generate_hmm_model(model, 1)

/*
init_set(double ***set, hmm_seq *seq, hmm_model *model);
void init_set_v(double **set_p, hmm_seq *seq);
void init_ksiset(double ****ksiset, hmm_seq *seq, hmm_model *model);
void init_gammaset(double ***gammaset_p, hmm_seq *seq, hmm_model *model);

void free_set_v(double *set);
void free_set(double **set,  hmm_seq *seq);
void free_ksiset(double ***ksiset, hmm_seq *seq, hmm_model *model);
void free_gammaset(double **gammaset, hmm_seq *seq);
*/

/**
 *
 * @param seq
 * @param t
 * @param n
 * @param m
 * @param type
 * @return
 */
byte generate_hmm_seq(hmm_seq *seq, hmm_model *model);
byte generate_random_hmm_seq(hmm_seq *seq, byte m);



byte forward_algorithm(hmm_seq *seq, hmm_model *model, double **set, double *set_v);
byte backward_algorithm(hmm_seq *seq, hmm_model *model, double **set, double *set_v);
double estimation_sequence_forward(hmm_seq *seq, hmm_model *model, double **set, double *set_v);
void double_probability_norm(hmm_seq *seq, hmm_model *model, double estimation_seq, double **alphaset,
                             double *alhaset_v, double **betaset, double *betaset_v, double  ***ksiset);
void marginaol_probability_norm(hmm_seq *seq, hmm_model *model, double estimation_seq, double **alphaset,
                                double *alhaset_v, double **betaset, double *betaset_v, double  **gammaset);

void estimation_model(hmm_seq *seq, hmm_model *model, double eps, double *likehood);
void estimation_model_gl(hmm_seq *seq, hmm_model *model, uint32_t iter, double eps); //search for global maximum of likehood.

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_HMM_H
