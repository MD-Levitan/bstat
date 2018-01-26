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

typedef struct sequence_cmm_s{
    cmm_model *model;   /**< CMM */
    qword T;            /**< The length of sequence. */
    byte  m;            /**< The size of alphabet of seq, if model == NULL. */
    byte* array;        /**< The pointer to sequence. */
} cmm_seq;

void init_cmm_model(cmm_model *ctx, byte n);
void init_cmm_seq(cmm_seq *ctx, qword t, cmm_model *ctx2);
void copy_cmm_model(cmm_model *dist, cmm_model *src);
void free_cmm_model(cmm_model *ctx);
void free_cmm_seq(cmm_seq *ctx);

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

byte generate_cmm_seq(cmm_seq *seq);
byte generate_random_cmm_seq(cmm_seq *seq, byte m);


byte MLE_algorithm(cmm_seq *seq, cmm_model *model);
byte bootstrap(cmm_seq *seq, cmm_model *model, word repeats);
byte smoothed_estimators(cmm_seq *seq, cmm_model *model, word repeats, double u);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_CMM_H
