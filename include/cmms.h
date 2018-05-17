#ifndef __BSTAT_CMM_S_H
#define __BSTAT_CMM_S_H

#include "bstat.h"

#ifdef __cplusplus
extern "C"{
#endif

/**
 * @brief Params of CMM_S(Chain Markov Model).
 *
 *
 */
typedef struct cmms{
    byte S;             /**< Order of CMM(size of memory). */
    byte N;             /**< Size of the set of states. */
    double *Pi;         /**< Initial vector. */
    double **P;         /**< One-step transition matrix. */
} cmms_model;


void init_cmms_model(cmms_model *ctx, byte n, byte s);
void copy_cmms_model(cmms_model *dist, cmms_model *src);
void set_cmms_model(cmms_model *dist, double *Pi, double **P);
void free_cmms_model(cmms_model *ctx);

/**
 * @brief Generate params of CMM_S.
 *
 * @param model The pointer of model.
 * @param type  The type of generation params of CMM_S.
 *              @p type = {0, 1}, where 0 - random, 1 - uniform.
 * @return      Return #Success or #Error.
 *
 */
byte generate_cmms_model(cmms_model *model, byte type);
#define generate_random_cmms_model(model) generate_cmms_model(model, 0)
#define generate_uniform_cmms_model(model) generate_cmms_model(model, 1)

/**
 * Generation of sequence.
 * 1. Using CMMs model.
 */
byte generate_cmms_sequence(sequence *seq, cmms_model *model);

/**
 * Algorithms of estimation.
 * Type: 1. MLE; 2. Bootstrap; 3. Smoothed alg.
 * Type of data: 1. sequence; 2. stream;
 */
byte MLE_algorithm_s(sequence *seq, cmms_model *model, double *n);
byte bootstrap_s(sequence *seq, cmms_model *model, word repeats, double *n);
byte smoothed_estimators_s(sequence *seq, cmms_model *model, word repeats, double u, double *n);

byte MLE_algorithm_s_stream(istream *str, cmms_model *model, double *n);
byte bootstrap_s_stream(istream *str, cmms_model *model, word repeats, double *n);
byte smoothed_estimators_s_stream(istream *str, cmms_model *model, word repeats, double u, double *n);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_CMM_S_H
