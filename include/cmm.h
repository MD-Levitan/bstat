#ifndef __BSTAT_CMM_H
#define __BSTAT_CMM_H

#include "bstat.h"

#ifdef __cplusplus
extern "C"{
#endif

/**
 * @brief Params of CMM(Chain Markov Model).
 *
 *
 */
typedef struct cmm_s{
    byte N;             /**< Size of the set of states. */
    double *Pi;         /**< Initial vector. */
    double **P;         /**< One-step transition matrix. */
} cmm_model;


void init_cmm_model(cmm_model *ctx, byte n);
void copy_cmm_model(cmm_model *dist, cmm_model *src);
void free_cmm_model(cmm_model *ctx);

/**
 * @brief Generate params of CMM.
 *
 * @param model The pointer of model.
 * @param type  The type of generation params of CMM.
 *              @p type = {0, 1}, where 0 - random, 1 - uniform.
 * @return      Return #Success or #Error.
 *
 */
byte generate_cmm_model(cmm_model *model, byte type);
#define generate_random_cmm_model(model) generate_cmm_model(model, 0)
#define generate_uniform_cmm_model(model) generate_cmm_model(model, 1)

byte generate_cmm_sequence(sequence *seq, cmm_model *model);
byte generate_random_sequence(sequence *seq, byte m);



byte MLE_algorithm_stream(stream *str, cmm_model *model, double *n);
byte MLE_algorithm(sequence *seq, cmm_model *model, double *n);
byte bootstrap(sequence *seq, cmm_model *model, word repeats, double *n);
byte smoothed_estimators(sequence *seq, cmm_model *model, word repeats, double u, double *n);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_CMM_H
