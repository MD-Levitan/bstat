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

byte generate_hmm_model(hmm_model *model,  uint8_t type, uint32_t* param){
    if(_entropy == NULL || model == NULL)
        return ERROR;

    if(type != 0 && type !=1 && type != 2 )
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
        /*
         *  index1 - index in surface of size(N * (N-1)), which represent matrix P.
         *  index2 - index in surface of size(N * (M-1)), which represent matrix C.
         *  num1 - number of pieces of surface(P).
         *  num2 - number of pieces of surface(C).
         *
         *  1/6 * N * (N + 1) * (N + 2) = sum 1/2 * (n^2 + n) from 1 to N       --      (formula 1).
         *
         */


        uint32_t num1 = param[0];
        uint32_t num2 = param[1];

        double step1 = 1.0 / (num1 - 1);
        double step2 = 1.0 / (num2 - 1);

        uint64_t index1 = param[2];
        uint64_t index2 = param[3];

        //check if indexes are in correct range.
        if(index1 > pow(num1, model->N * (model->N - 1)) || index2 > pow(num1, model->N * (model->M - 1))){
            return ERROR;
        }
        
        //Set PI
        for (byte i = 0; i < model->N; ++i) {
            model->Pi[i] = (double) 1 / model->N;
        }

        //set P
        {
            //calculate using formula 1 for matrix P.
            double max_size_double = 1;
            for (uint8_t i = 1; i <= model->N - 1; ++i) {
                max_size_double *= ((num1 - 1) + i) / i;
            }
            uint32_t max_size = (uint32_t) max_size_double;
            uint32_t num_m1 = pow(num1, model->N - 1);
        
            for (byte i = 0; i < model->N; ++i) {
                uint32_t v = (index1 % num_m1);
                index1 /= num_m1;
                double res = 1;
                uint32_t current_size = 0;
                for (byte j = 0; j < model->N - 1; ++j) {
                    current_size += v % num1;
                    if (current_size >= max_size)
                        model->P[i][j] = 0;
                    double value = step1 * (v % num1);
                    res -= value;
                    if (res < 0)
                        return ERROR;
                    model->P[i][j] = value;
                    v /= num1;
                }
                model->P[i][model->N - 1] = res;
            }
        }

        //set C
        {
            //calculate using formula 1 for matrix P.
            double max_size_double = 1;
            for (uint8_t i = 1; i <= model->M - 1; ++i) {
                max_size_double *= (num2 + i) / i;
            }
            uint32_t max_size = (uint32_t) max_size_double;
        
            for (byte i = 0; i < model->N; ++i) {
                uint32_t v = index2 % num1;
                index2 /= num2;
                double res = 1;
                uint32_t current_size = 0;
                for (byte j = 0; j < model->M - 1; ++j) {
                    current_size += v % num2;
                    if (current_size >= max_size)
                        model->C[i][j] = 0;
                    double value = step2 * (v % num2);
                    res -= value;
                    if (res < 0)
                        return ERROR;
                    model->C[i][j] = value;
                    v /= num1;
                }
                model->C[i][model->M - 1] = res;
            }
        }
    }

    return SUCCESS;
}

byte generate_hmm_sequence(sequence *seq,  hmm_model *model) {
    if(_entropy == NULL || seq == NULL || model == NULL)
        return ERROR;

    byte hidden = ch(model->Pi, model->N);
    byte visible = ch(model->C[hidden], model->M);
    for (qword i = 0; i < seq->T; ++i) {
        seq->array[i] = visible;
        hidden = ch(model->P[hidden], model->N);
        visible = ch(model->C[hidden], model->M);
    }
}


//////////////////////////////////////////////////////////////////////
/////////// Free and init functions for intermediate vars. ///////////
//////////////////////////////////////////////////////////////////////

/*
 * set = {alphaset, betaset};
 * set_v = {alphaset_v, betaset_v};
 */

void init_set(double ***set_p, sequence *seq, hmm_model *model){
    if(!seq || !model)
        return;
    double **set;
    _memcheck(set, seq->T * sizeof(double*));
    for (qword t = 0; t < seq->T; ++t) {
        _memcheck(set[t], model->N * sizeof(double));
    }
    *set_p = set;
}

void free_set(double **set,  sequence *seq) {
    if (!set || !seq)
        return;

    for (qword t = 0; t < seq->T; ++t) {
        free(set[t]);
    }
    free(set);
}

void init_set_v(double **set_p, sequence *seq) {
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

void init_ksiset(double ****ksiset_p, sequence *seq, hmm_model *model){
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

void free_ksiset(double ***ksiset, sequence *seq, hmm_model *model) {
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

void init_gammaset(double ***gammaset_p, sequence *seq, hmm_model *model) {
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

void free_gammaset(double **gammaset, sequence *seq){
    if(!gammaset || !seq)
        return;
    for (qword i = 0; i < seq->T; ++i)
        free(gammaset[i]);
    free(gammaset);
}
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
byte forward_algorithm_(sequence *seq, hmm_model *model, double **set) {
    if(seq == NULL || model == NULL || seq->m != model->M)
        return ERROR;
    double **alpha_set = set;

    for (byte i = 0; i < model->N; ++i) {
        alpha_set[0][i] = model->C[i][seq->array[0]] * model->Pi[i];
    }

    for (qword t = 1; t < seq->T; ++t) {
        for (byte i = 0; i < model->N; ++i) {
            double sum = 0;
            for (byte j = 0; j < model->N; ++j) {
                sum += model->P[j][i] * alpha_set[t - 1][j];
            }
            alpha_set[t][i] = model->C[i][seq->array[t]] * sum;
        }
    }
    return SUCCESS;

}

byte forward_algorithm(sequence *seq, hmm_model *model, double **set, double *set_v){
    if(seq == NULL || model == NULL || seq->m != model->M)
        return ERROR;
    double *alpha, alpha_v = 0;
    double **alpha_set = set, *alpha_v_set = set_v;
    _memcheck(alpha, model->N * sizeof(double));
    

    for (byte i = 0; i < model->N; ++i) {
        alpha[i] = model->C[i][seq->array[0]] * model->Pi[i];
        if(isnan(alpha[i])){
            printf("adsd");
        }
        alpha_v += alpha[i];
    }
    alpha_v /= model->N;
    alpha_v_set[0] = alpha_v;
    alpha_v = !alpha_v ? 1 : alpha_v; //if alpha_v = 0;
    for (byte i = 0; i < model->N; ++i) {
        alpha[i] /= alpha_v;
        if(isnan(alpha[i])){
            printf("adsd");
        }
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
        alpha_v = !alpha_v ? 1 : alpha_v; //if alpha_v = 0;
        for (byte i = 0; i < model->N; ++i) {
            alpha[i] /= alpha_v;
            alpha_set[t][i] = alpha[i];
            if(isnan(alpha_set[t][i])){
                printf("adsd");
            }
        }
    }
    free(alpha);
    return SUCCESS;
}

byte backward_algorithm(sequence *seq, hmm_model *model, double **set, double *set_v){
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
        beta_v  /= model->N;
        beta_v_set[t - 1] = beta_v;
        beta_v = !beta_v ? 1 : beta_v;
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
        beta_v = !beta_v ? 1 : beta_v;
        for (byte i = 0; i < model->N; ++i) {
            beta_set[0][i] /= beta_v;
        }
    }

    return SUCCESS;
}

byte backward_algorithm_(sequence *seq, hmm_model *model, double **set){
    if(seq == NULL || model == NULL || seq->m != model->M)
        return ERROR;
    double **beta_set = set;

    for (byte i = 0; i < model->N; ++i)
        beta_set[seq->T - 1][i] = 1;

    for (qword t = seq->T - 1; t >= 2; --t) {
        for (byte i = 0; i < model->N; ++i) {
            double sum = 0;
            for (byte j = 0; j < model->N; ++j) {
                sum += model->P[i][j] * model->C[j][seq->array[t]] * beta_set[t][j];
            }
            beta_set[t - 1][i] = sum;
        }
    }

    //t = 1
    {
        for (byte i = 0; i < model->N; ++i) {
            double sum = 0;
            for (byte j = 0; j < model->N; ++j) {
                sum += model->P[i][j] * model->C[j][seq->array[1]] * beta_set[1][j];
            }
            beta_set[0][i] = sum;
        }
    }

    return SUCCESS;
}

double estimation_sequence_forward(sequence *seq, hmm_model *model, double **set, double *set_v){
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

double estimation_sequence_forward_(sequence *seq, hmm_model *model, double **set){
    double est = 0;
    for (byte i = 0; i < model->N; ++i) {
        est += set[seq->T - 1][i];
    }
    return est;
}

double estimation_sequence_forward_backward(sequence *seq, double *set_v, double **set){
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

void double_probability(sequence *seq, hmm_model *model, double estimation_seq, double **alphaset, double *alhaset_v,
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

void marginal_probability(sequence *seq, hmm_model *model, double estimation_seq, double **alphaset, double *alhaset_v,
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
void estimation_model(sequence *seq, hmm_model *model_i, double eps, double *likehood) {
    
    ///Init block
    hmm_model model;
    double est, prev_est;
    //double std_deviation = 0, prev_std_deviation = 1;
    double **alphaset;
    double **betaset;
    double *alphaset_v;
    double *betaset_v;
    double ***ksiset;
    double **gammaset;
    
    init_hmm_model(&model, model_i->N, model_i->M);
    copy_hmm_model(&model, model_i);
    //generate_random_hmm_model(&model);
    init_set(&alphaset, seq, model_i);
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
    if (isnan(prev_est) || isinf(prev_est)) {
        if (likehood != NULL)
            *likehood = prev_est;
        copy_hmm_model(model_i, &model);
        free_hmm_model(&model);
        return;
    }
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
                if (isnan(new_model.P[i][j]))
                    new_model.P[i][j] = 0;
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
                if (isnan(new_model.C[i][j]))
                    new_model.C[i][j] = 0;
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
        
        if (isnan(prev_est) || isnan(est)) {
            printf("dfsdf");
        }
        if (fabs(prev_est - est) <= eps)
            break;
        
        prev_est = est;
    }
    if (likehood != NULL)
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

void estimation_model_gl(sequence *seq, hmm_model *model_i, uint32_t iter, double eps, double *likehood) {
    double current_est, best_est = current_est = 0;
    hmm_model model;
    init_hmm_model(&model, model_i->N, model_i->M);
    uint64_t parts1 = pow(iter, model_i->N * (model_i->N - 1));
    uint64_t parts2 = pow(iter, model_i->N * (model_i->M - 1));
    double est = log(0);
    for (uint64_t i = 0; i < parts1; ++i) {
        for (uint64_t j = 0; j < parts2; ++j) {
            hmm_model model_iter;
            init_hmm_model(&model_iter, model_i->N, model_i->M);
            uint32_t params[4] = {iter, iter, i, j};
            if(generate_hmm_model(&model_iter, 2, params) == ERROR)
                continue;

//            printf("NEW MODEL\n");
//            printf("Pi:  ");
//            for (int j = 0; j < model_iter.N; ++j) {
//                printf("%lf ", model_iter.Pi[j]);
//            }
//
//            printf("\nP:  ");
//            for (byte i = 0; i < model_iter.N; ++i) {
//                for (int j = 0; j < model_iter.N; ++j) {
//                    printf("%lf ", model_iter.P[i][j]);
//                }
//                printf("\n\t");
//            }
//
//            printf("\nC:  ");
//            for (byte i = 0; i < model_iter.N; ++i) {
//                for (int j = 0; j < model_iter.M; ++j) {
//                    printf("%lf ", model_iter.C[i][j]);
//                }
//                printf("\n\t");
//            }
    
    
            double current_est = 0;
            estimation_model(seq, &model_iter, 0.01, &current_est);
           // printf("estimation = %f", current_est);
            if(est < current_est){
                copy_hmm_model(&model, &model_iter);
                est = current_est;
            }
         //   printf("\n\n====================================================\n");
        }
    }
    if(likehood != NULL)
        *likehood = est;
    printf("NEW GLOBAL EST = %f \n", est);
    copy_hmm_model(model_i, &model);
    free_hmm_model(&model);
}