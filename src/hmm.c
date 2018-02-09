#include <hmm.h>
#include <libnet.h>
#include "hmm.h"

void init_hmm_model(hmm_model *ctx, byte n, byte m){
    if(!ctx)
        return;
    ctx->M = m;
    ctx->N = n;

    _memcheck(ctx->P,  ctx->N * sizeof(double*));
    _memcheck(ctx->Pi, ctx->N * sizeof(double));
    _memcheck(ctx->C,  ctx->N * sizeof(double*));

    for (byte i = 0; i < ctx->N; ++i) {
        _memcheck(ctx->P[i], ctx->N * sizeof(double));
        _memcheck(ctx->C[i], ctx->M * sizeof(double));
    }
}

void init_hmm_seq(hmm_seq *ctx, qword t, hmm_model *ctx2){
    if(!ctx)
        return;

    ctx->T = t;
    _memcheck(ctx->array, ctx->T * sizeof(byte));
    if(ctx2 != NULL) {
        ctx->m = ctx2->M;
        ctx->model = ctx2;
    }
}

void copy_hmm_model(hmm_model *dist, hmm_model *src){
    if(!dist || !src || dist->M != src->M || dist->N != src->N)
        return;

    for (int j = 0; j < src->N; ++j) {
        dist->Pi[j] = src->Pi[j];
    }

    for (byte i = 0; i < src->N; ++i) {
        for (int j = 0; j < src->N; ++j) {
            dist->P[i][j] = src->P[i][j];
        }
    }

    for (byte i = 0; i < src->N; ++i) {
        for (int j = 0; j < src->M; ++j) {
            dist->C[i][j] = src->C[i][j];
        }
    }
}


void free_hmm_model(hmm_model *ctx){
    if(!ctx)
        return;

    for (byte i = 0; i < ctx->N; ++i) {
        free(ctx->P[i]);
        free(ctx->C[i]);
    }

    free(ctx->Pi);
    free(ctx->P);
    free(ctx->C);

    ctx->N = 0;
    ctx->M = 0;
}

void free_hmm_seq(hmm_seq *ctx){
    if(!ctx)
        return;

    free(ctx->array);
    ctx->T = 0;
}

byte generate_hmm_model(hmm_model *model,  uint8_t type, uint32_t* param){
    if(_entropy == NULL || model == NULL)
        return ERROR;

    if(type != 0 && type !=1 && type != 2)
        return ERROR;

    if(type == 0) {
        //Set PI
        {
            byte *arr;
            _memcheck(arr, model->N * sizeof(byte));
            qword sum = entropy_s(arr, model->N);
            for (byte i = 0; i < model->N; ++i) {
                model->Pi[i] = (double) arr[i] / sum;
            }
            free(arr);
        }

        //set P
        for (byte i = 0; i < model->N; ++i) {

            byte *arr;
            _memcheck(arr, model->N * sizeof(byte));
            qword sum = entropy_s(arr, model->N);
            for (byte j = 0; j < model->N; ++j)
                model->P[i][j] = (double) arr[j] / sum;

            free(arr);
        }

        //set C
        for (byte i = 0; i < model->N; ++i) {

            byte *arr;

            _memcheck(arr, model->M * sizeof(byte));
            qword sum = entropy_s(arr, model->M);
            for (byte j = 0; j < model->M; ++j)
                model->C[i][j] = (double) arr[j] / sum;

            free(arr);
        }
    }else if(type == 1){
        //Set PI
        for (byte i = 0; i < model->N; ++i) {
            model->Pi[i] = (double) 1 / model->N;
        }

        //set P
        for (byte i = 0; i < model->N; ++i) {
            for (byte j = 0; j < model->N; ++j)
                model->P[i][j] = (double) 1 / model->N;
        }

        //set C
        for (byte i = 0; i < model->N; ++i) {
            for (byte j = 0; j < model->M; ++j)
                model->C[i][j] = (double) 1 / model->M;
        }
    }else{
        if(param == NULL)
            return ERROR;
        uint32_t num1 = param[0];
        uint32_t num2 = param[1];
        uint32_t index = param[2];
        uint32_t index2 = index / num1;
        uint32_t index1 = index % num1;
        

        //Set PI
        for (byte i = 0; i < model->N; ++i) {
            model->Pi[i] = (double) 1 / model->N;
        }

        //set P
        for (byte i = 0; i < model->N; ++i) {
            for (byte j = 0; j < model->N; ++j)
                model->P[i][j] = (double) index1 / model->N;
        }

        //set C
        for (byte i = 0; i < model->N; ++i) {
            for (byte j = 0; j < model->M; ++j)
                model->C[i][j] = (double) 1 / model->M;
        }
    }

    return SUCCESS;
}

byte generate_random_hmm_seq(hmm_seq *seq, byte m){
    if(_entropy == NULL || seq == NULL)
        return ERROR;

    seq->m = m;
    entropy_m(seq->array, seq->T, seq->m);
}

byte generate_hmm_seq(hmm_seq *seq,  hmm_model *model) {
    if(_entropy == NULL || seq == NULL || model == NULL)
        return ERROR;

    byte hidden = ch(seq->model->Pi, seq->model->N);
    byte visible = ch(seq->model->C[hidden], seq->model->M);
    for (qword i = 0; i < seq->T; ++i) {
        seq->array[i] = visible;
        hidden = ch(seq->model->P[hidden], seq->model->N);
        visible = ch(seq->model->C[hidden], seq->model->M);
    }
}


//////////////////////////////////////////////////////////////////////
/////////// Free and init functions for intermediate vars. ///////////
//////////////////////////////////////////////////////////////////////

/*
 * set = {alphaset, betaset};
 * set_v = {alphaset_v, betaset_v};
 */

void init_set(double ***set_p, hmm_seq *seq, hmm_model *model){
    if(!seq || !model)
        return;
    double **set;
    _memcheck(set, seq->T * sizeof(double*));
    for (qword t = 0; t < seq->T; ++t) {
        _memcheck(set[t], model->N * sizeof(double));
    }
    *set_p = set;
}

void free_set(double **set,  hmm_seq *seq) {
    if (!set || !seq)
        return;

    for (qword t = 0; t < seq->T; ++t) {
        free(set[t]);
    }
    free(set);
}

void init_set_v(double **set_p, hmm_seq *seq) {
    if (!seq)
        return;
    double *set;
    _memcheck(set, seq->T * sizeof(double));
    *set_p = set;
}

void free_set_v(double *set){
    if(!set)
        return;
    free(set);
}

void init_ksiset(double ****ksiset_p, hmm_seq *seq, hmm_model *model){
    if(!seq || !model)
        return;
    double ***ksiset;
    _memcheck(ksiset, (seq->T - 1) * sizeof(double**));
    for(qword i = 0; i < seq->T - 1; ++i) {
        _memcheck(ksiset[i], model->N * sizeof(double *));
        for(byte k = 0; k < model->N; ++k){
            _memcheck(ksiset[i][k], model->N * sizeof(double));
        //    memset(ksiset[i][k], 0, sizeof(double) * model->N);//not good solution
        }
    }
    *ksiset_p = ksiset;
}

void free_ksiset(double ***ksiset, hmm_seq *seq, hmm_model *model) {
    if(!ksiset || !seq || !model)
        return;

    for(qword i = 0; i < seq->T - 1; ++i) {
        for(byte k = 0; k < model->N; ++k){
            free(ksiset[i][k]);
        }
        free(ksiset[i]);
    }
    free(ksiset);
}

void init_gammaset(double ***gammaset_p, hmm_seq *seq, hmm_model *model) {
    if(!seq || !model)
        return;
    double **gammaset;
    _memcheck(gammaset, seq->T * sizeof(double *));
    for (qword i = 0; i < seq->T; ++i) {
        _memcheck(gammaset[i], model->N * sizeof(double));
        //memset(gammaset[i], 0, sizeof(double) * model->N);//not good solution
    }
    *gammaset_p = gammaset;
}

void free_gammaset(double **gammaset, hmm_seq *seq){
    if(!gammaset || !seq)
        return;
    for (qword i = 0; i < seq->T; ++i)
        free(gammaset[i]);
    free(gammaset);
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

byte forward_algorithm(hmm_seq *seq, hmm_model *model, double **set, double *set_v){
    if(seq == NULL || model == NULL || seq->m != model->M)
        return ERROR;
    double *alpha, alpha_v = 0;
    double **alpha_set = set, *alpha_v_set = set_v;
    _memcheck(alpha, model->N * sizeof(double));

    for (byte i = 0; i < model->N; ++i) {
        alpha[i] = model->C[i][seq->array[0]] * model->Pi[i];
        alpha_v += alpha[i];
    }
    alpha_v /= model->N;
    alpha_v_set[0] = alpha_v;
    for (byte i = 0; i < model->N; ++i) {
        alpha[i] /= alpha_v;
        alpha_set[0][i] = alpha[i];
    }

    for (qword t = 1; t < seq->T; ++t) {
        alpha_v = 0;
        for (byte i = 0; i < model->N; ++i) {
            double sum = 0;
            for (byte j = 0; j < model->N; ++j) {
                sum += model->P[j][i] * alpha_set[t - 1][j];
            }
            alpha[i] = model->C[i][seq->array[t]] * sum;
            alpha_v += alpha[i];
        }
        alpha_v /= model->N;
        alpha_v_set[t] = alpha_v;
        for (byte i = 0; i < model->N; ++i) {
            alpha[i] /= alpha_v;
            alpha_set[t][i] = alpha[i];
        }
    }
    free(alpha);
    return SUCCESS;
}

byte backward_algorithm(hmm_seq *seq, hmm_model *model, double **set, double *set_v){
    if(seq == NULL || model == NULL || seq->m != model->M)
        return ERROR;
    double beta_v = 1;
    double **beta_set = set, *beta_v_set = set_v;


    beta_v_set[seq->T - 1] = beta_v;
    for (byte i = 0; i < model->N; ++i)
        beta_set[seq->T - 1][i] = 1;

    for (qword t = seq->T - 1; t >= 2; --t) {
        beta_v = 0;
        for (byte i = 0; i < model->N; ++i) {
            double sum = 0;
            for (byte j = 0; j < model->N; ++j) {
                sum += model->P[i][j] * model->C[j][seq->array[t]] * beta_set[t][j];
            }
            beta_set[t - 1][i] = sum;
            beta_v += beta_set[t - 1][i];
        }
        beta_v /= model->N;
        beta_v_set[t - 1] = beta_v;
        for (byte i = 0; i < model->N; ++i) {
            beta_set[t - 1][i] /= beta_v;
        }
    }

    //t = 1
    {
        beta_v = 0;
        for (byte i = 0; i < model->N; ++i) {
            double sum = 0;
            for (byte j = 0; j < model->N; ++j) {
                sum += model->P[i][j] * model->C[j][seq->array[1]] * beta_set[1][j];
            }
            beta_set[0][i] = sum;
            beta_v += beta_set[0][i];
        }
        beta_v /= model->N;
        beta_v_set[0] = beta_v;
        for (byte i = 0; i < model->N; ++i) {
            beta_set[0][i] /= beta_v;
        }
    }

    return SUCCESS;
}

double estimation_sequence_forward(hmm_seq *seq, hmm_model *model, double **set, double *set_v){
    double est = 0;
    for (byte i = 0; i < model->N; ++i) {
        est += set[seq->T - 1][i];
    }
    est = log2(est);
    for (qword j = 0; j < seq->T; ++j) {
        est += log2(set_v[j]);
    }
    return est;
}

double estimation_sequence_forward_backward(hmm_seq *seq, double *set_v, double **set){
    /*double est = 0;
    for (byte i = 0; i < seq->m; ++i) {
        est += set[seq->T - 1][i];
    }
    est = log2(est);
    for (qword j = 0; j < seq->T; ++j) {
        est += log2(set_v[j]);
    }
    return est;*/
}

void double_probability(hmm_seq *seq, hmm_model *model, double estimation_seq, double **alphaset, double *alhaset_v,
                             double **betaset, double *betaset_v, double  ***ksiset) {

    for (qword t = 0; t < seq->T - 1; ++t) {
        for (byte i = 0; i < model->N; ++i) {
            for (byte j = 0; j < model->N; ++j) {
                double sum = log2(alphaset[t][i]);
                for (qword teta = 0; teta <= t; ++teta)
                    sum += log2(alhaset_v[teta]);
                sum += log2(model->P[i][j]);
                sum += log2(model->C[j][seq->array[t + 1]]);
                sum += log2(betaset[t+1][j]);
                for (qword teta = t; teta < seq->T - 1; ++teta)
                    sum += log2(betaset_v[teta + 1]);
                sum -= estimation_seq;
                ksiset[t][i][j] = exp2(sum);
            }
        }
    }
}

void marginal_probability(hmm_seq *seq, hmm_model *model, double estimation_seq, double **alphaset, double *alhaset_v,
                          double **betaset, double *betaset_v, double **gammaset){
    for(qword t = 0; t < seq->T - 1; ++t) {
        for(byte i = 0; i < model->N; ++i) {
            double sum = log2(alphaset[t][i]);
            for (qword teta = 0; teta <= t; ++teta)
                sum += log2(alhaset_v[teta]);
            sum += log2(betaset[t][i]);
            for (qword teta = t; teta < seq->T; ++teta)
                sum += log2(betaset_v[teta]);
            sum -= estimation_seq;
            gammaset[t][i] = exp2(sum);
        }
    }
}


/*
 * Need to repair pointer.
 * Try to add search of global maximum.
 */
void estimation_model(hmm_seq *seq, hmm_model *model_i, double eps, double *likehood){

    ///Init block
    hmm_model model;
    double est, prev_est;
    //double std_deviation = 0, prev_std_deviation = 1;
    double ** alphaset;
    double ** alphaset_new;
    double ** betaset;
    double * alphaset_v;
    double * betaset_v;
    double *** ksiset;
    double **gammaset;

    init_hmm_model(&model, model_i->N, model_i->M);
    generate_random_hmm_model(&model);
    init_set(&alphaset, seq, model_i);
    init_set(&alphaset_new, seq, model_i);
    init_set(&betaset, seq, model_i);
    init_set_v(&alphaset_v, seq);
    init_set_v(&betaset_v, seq);
    init_ksiset(&ksiset, seq, model_i);
    init_gammaset(&gammaset, seq, model_i);
    /////

    /// First iteration
    forward_algorithm(seq, &model, alphaset, alphaset_v);
    backward_algorithm(seq, &model, betaset, betaset_v);
    prev_est = estimation_sequence_forward(seq, &model, alphaset, alphaset_v);
    ///

    while (1) {
        hmm_model new_model;
        init_hmm_model(&new_model, model.N, model.M);

        marginal_probability(seq, &model, prev_est, alphaset, alphaset_v, betaset, betaset_v, gammaset);
        double_probability(seq, &model, prev_est, alphaset, alphaset_v, betaset, betaset_v, ksiset);

        for (byte i = 0; i < model.N; ++i)
            new_model.Pi[i] = gammaset[0][i];

        for (byte i = 0; i < model.N; ++i)
            for (byte j = 0; j < model.N; ++j) {
                double sum1, sum2 = sum1 = 0;
                for (qword t = 0; t < seq->T - 1; ++t) {
                    sum1 += ksiset[t][i][j];
                    sum2 += gammaset[t][i];
                }
                new_model.P[i][j] = sum1 / sum2;
            }

        for (byte i = 0; i < model.N; ++i)
            for (byte j = 0; j < seq->m; ++j) {
                double sum1, sum2 = sum1 = 0;
                for (qword t = 0; t < seq->T - 1; ++t) {
                    if (seq->array[t] == j)
                        sum1 += gammaset[t][i];
                    sum2 += gammaset[t][i];
                }
                new_model.C[i][j] = sum1 / sum2;
            }

        /// Calculate likehood for new_model
        forward_algorithm(seq, &new_model, alphaset, alphaset_v);
        backward_algorithm(seq, &new_model, betaset, betaset_v);
        est = estimation_sequence_forward(seq, &new_model, alphaset, alphaset_v);
        ///

        // std_deviation += standart_deviation_matrix(model.P, new_model.P, model.N, model.N);
        // std_deviation += standart_deviation_matrix(model.C, new_model.C, model.N, model.M);

        copy_hmm_model(&model, &new_model);
        free_hmm_model(&new_model);

        if (est < prev_est)
            fprintf(stderr, "%s. Likehood of prev. iter.(%f) < Likehood of new iter.(%f).\n",
                    "Warning! Algorithm EM doesn't work properly!", prev_est, est);

        if (fabs(prev_est - est) <= eps)
            break;

        prev_est = est;
    }
    if(likehood != NULL)
        *likehood = est;
    copy_hmm_model(model_i, &model);
    
    ///free block
    free_set(alphaset, seq);
    free_set(betaset, seq);
    free_set_v(alphaset_v);
    free_set_v(betaset_v);
    free_ksiset(ksiset, seq, &model);
    free_gammaset(gammaset, seq);
    free_hmm_model(&model);
    ////
}

void estimation_model_gl(hmm_seq *seq, hmm_model *model_i, uint32_t iter, double eps) {
    double current_est, best_est = current_est = 0;
    hmm_model model;
    init_hmm_model(&model, model_i->N, model_i->M);


    copy_hmm_model(model_i, &model);
    free_hmm_model(&model);
}