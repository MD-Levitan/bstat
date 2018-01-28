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

typedef struct sequence_cmms_s{
    cmms_model *model;   /**< CMM */
    qword T;            /**< The length of sequence. */
    byte  m;            /**< The size of alphabet of seq, if model == NULL. */
    byte* array;        /**< The pointer to sequence. */
} cmms_seq;

void init_cmms_model(cmms_model *ctx, byte n, byte s);
void init_cmms_seq(cmms_seq *ctx, qword t, cmms_model *ctx2);
void copy_cmms_model(cmms_model *dist, cmms_model *src);
void free_cmms_model(cmms_model *ctx);
void free_cmms_seq(cmms_seq *ctx);

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

byte generate_cmms_seq(cmms_seq *seq);
byte generate_random_cmms_seq(cmms_seq *seq, byte m);


byte MLE_algorithm(cmms_seq *seq, cmms_model *model, double *n);
byte bootstrap(cmms_seq *seq, cmms_model *model, word repeats, double *n);
byte smoothed_estimators(cmms_seq *seq, cmms_model *model, word repeats, double u, double *n);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif //__BSTAT_CMM_S_H
