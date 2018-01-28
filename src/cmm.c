#include <libnet.h>
#include "cmm.h"

/////////////////////////////////////////////////
//////////////// INIT, FREE, COPY ///////////////
/////////////////////////////////////////////////

void init_cmm_model(cmm_model *ctx, byte n){
    if(!ctx)
        return;
    ctx->N = n;

    _memcheck(ctx->P,  ctx->N * sizeof(double*));
    _memcheck(ctx->Pi, ctx->N * sizeof(double));

    for (byte i = 0; i < ctx->N; ++i)
        _memcheck(ctx->P[i], ctx->N * sizeof(double));
}

void init_cmm_seq(cmm_seq *ctx, qword t, cmm_model *ctx2){
    if(!ctx)
        return;

    ctx->T = t;
    _memcheck(ctx->array, ctx->T * sizeof(byte));
    if(ctx2 != NULL) {
        ctx->model = ctx2;
        ctx->m = ctx2->N;
    }
}

void copy_cmm_model(cmm_model *dist, cmm_model *src){
    if(!dist || !src || dist->N != src->N)
        return;

    for (int j = 0; j < src->N; ++j) {
        dist->Pi[j] = src->Pi[j];
    }

    for (byte i = 0; i < src->N; ++i) {
        for (int j = 0; j < src->N; ++j) {
            dist->P[i][j] = src->P[i][j];
        }
    }
}

void free_cmm_model(cmm_model *ctx){
    if(!ctx)
        return;

    for (byte i = 0; i < ctx->N; ++i)
        free(ctx->P[i]);

    free(ctx->Pi);
    free(ctx->P);

    ctx->N = 0;
}

void free_cmm_seq(cmm_seq *ctx){
    if(!ctx)
        return;

    free(ctx->array);
    ctx->T = 0;
}

/////////////////////////////////////////////////


byte generate_random_cmm_seq(cmm_seq *seq, byte m){
    if(_entropy == NULL || seq == NULL)
        return ERROR;

    seq->m = m;
    entropy_m(seq->array, seq->T, seq->m);
}

byte generate_cmm_seq(cmm_seq *seq){
    if(_entropy == NULL || seq == NULL || seq->model == NULL)
        return ERROR;

    byte visible = ch(seq->model->Pi, seq->model->N);
    for (qword i = 0; i < seq->T; ++i) {
        seq->array[i] = visible;
        visible = ch(seq->model->P[visible], seq->model->N);
    }
}

byte generate_cmm_model(cmm_model *model, byte type){
    if(_entropy == NULL || model == NULL)
        return ERROR;

    if(type != 0 && type !=1)
        return ERROR;

    if(!type) {
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
            for (byte j = 0; j < model->N; ++j) {
                model->P[i][j] = (double) arr[j] / sum;
            }
            free(arr);
        }
    }else{
        //Set PI
        for (byte i = 0; i < model->N; ++i) {
            model->Pi[i] = (double) 1 / model->N;
        }

        //set P
        for (byte i = 0; i < model->N; ++i) {
            for (byte j = 0; j < model->N; ++j)
                model->P[i][j] = (double) 1 / model->N;
        }
    }

    return SUCCESS;
}

byte MLE_algorithm(cmm_seq *seq, cmm_model *model, double *_n) {
    if (!model || !seq || (seq->m && seq->m != model->N))
        return ERROR;

    //calculate Pi
    double value = 1.0 / model->N;
    for (byte i = 0; i < model->N; ++i) {
        model->Pi[i] = value;
    }

    //calculate n_i
    //Complexity = O(N) + 0(T-1)
    double *n;
    _memcheck(n, model->N * sizeof(double));
    for (byte i = 0; i < model->N; ++i) {
        n[i] = 0;
    }

    for (qword t = 0; t < seq->T - 1; ++t) {
        n[seq->array[t]]++;
    }

    //calculate P
    //Complexity = O(N*N*(T-1))
//    for (byte i = 0; i < model->N; ++i) {
//        for (byte j = 0; j < model->N; ++j) {
//            value = 0;
//            for (qword t = 0; t < seq->T - 1; ++t) {
//                if((seq->array[t] == i) && (seq->array[t + 1] == j))
//                    value++;
//            }
//            model->P[i][j] = n[i] ? value / n[i] : 1.0 / model->N;
//        }
//    }

    //Complexity = O(2*N*N) + 0(T-1)
    for (byte i = 0; i < model->N; ++i) {
        for (byte j = 0; j < model->N; ++j) {
            model->P[i][j] = 0;
        }
    }
    for (qword t = 0; t < seq->T - 1; ++t) {
        model->P[seq->array[t]][seq->array[t + 1]] ++ ;
    }
    for (byte i = 0; i < model->N; ++i) {
        for (byte j = 0; j < model->N; ++j) {
            model->P[i][j] = n[i] ? model->P[i][j] / n[i] : 1.0 / model->N;
        }
    }
    if(_n)
        memcpy(_n, n, model->N * sizeof(double));
    free(n);
    return SUCCESS;
}

byte bootstrap(cmm_seq *seq, cmm_model *model, word repeats, double *_n){
    if (!model || !seq || (seq->m && seq->m != model->N))
        return ERROR;
    cmm_model iter_model;
    cmm_seq iter_seq;
    cmm_model average_model;
    cmm_model it_model;
    init_cmm_model(&iter_model, model->N);
    init_cmm_model(&it_model, model->N);
    init_cmm_model(&average_model, model->N);
    double *n, *n_average;
    _memcheck(n, model->N * sizeof(double));
    _memcheck(n_average, model->N * sizeof(double));


    MLE_algorithm(seq, &iter_model, n);
    memcpy(n_average, n, model->N * sizeof(double));
    copy_cmm_model(&average_model, &iter_model);
    init_cmm_seq(&iter_seq, seq->T, &iter_model);
    for (word it = 1; it < repeats; ++it) {
        generate_cmm_seq(&iter_seq);
        MLE_algorithm(&iter_seq, &it_model, n);
        for(byte i = 0; i < model->N; ++i){
            n_average[i] += n[i];
        }
        for(byte i = 0; i < model->N; ++i){
            for(byte j = 0; j < model->N; ++j){
                   average_model.P[i][j] += it_model.P[i][j];
            }
        }
    }
    for(byte i = 0; i < model->N; ++i){
        for(byte j = 0; j < model->N; ++j){
            average_model.P[i][j] /= repeats;
        }
    }
    for(byte i = 0; i < model->N; ++i){
        n_average[i] /= repeats;
    }
    if(_n)
        memcpy(_n, n, model->N * sizeof(double));
    copy_cmm_model(model, &average_model);
    free_cmm_model(&average_model);
    free_cmm_model(&iter_model);
    free_cmm_model(&it_model);
    free_cmm_seq(&iter_seq);
    return SUCCESS;
}

byte smoothed_estimators(cmm_seq *seq, cmm_model *model, word repeats, double u, double *_n) {
    if (!model || !seq || (seq->m && seq->m != model->N))
        return ERROR;
    cmm_model iter_model;
    cmm_seq iter_seq;
    cmm_model average_model;
    cmm_model it_model;
    init_cmm_model(&iter_model, model->N);
    init_cmm_model(&it_model, model->N);
    init_cmm_model(&average_model, model->N);
    double gamma = pow(seq->T,-u);
    double omega = 1 + gamma * model->N;
    double *n, *n_average;
    _memcheck(n, model->N * sizeof(double));
    _memcheck(n_average, model->N * sizeof(double));


    MLE_algorithm(seq, &iter_model, n);
    copy_cmm_model(&average_model, &iter_model);
    for(byte i = 0; i < model->N; ++i){
        for(byte j = 0; j < model->N; ++j){
            iter_model.P[i][j] = (iter_model.P[i][j] + gamma) / omega;
        }
    }
    init_cmm_seq(&iter_seq, seq->T, &iter_model);
    for (word it = 1; it < repeats; ++it) {
        generate_cmm_seq(&iter_seq);
        MLE_algorithm(&iter_seq, &it_model, n);
        for(byte i = 0; i < model->N; ++i){
            n_average[i] += n[i];
        }
        for(byte i = 0; i < model->N; ++i){
            for(byte j = 0; j < model->N; ++j){
                average_model.P[i][j] += it_model.P[i][j];
            }
        }
    }
    for(byte i = 0; i < model->N; ++i){
        for(byte j = 0; j < model->N; ++j){
            average_model.P[i][j] /= repeats;
        }
    }
    for(byte i = 0; i < model->N; ++i){
        n_average[i] /= repeats;
    }
    if(_n)
        memcpy(_n, n, model->N * sizeof(double));

    copy_cmm_model(model, &average_model);
    free_cmm_model(&average_model);
    free_cmm_model(&iter_model);
    free_cmm_model(&it_model);
    free_cmm_seq(&iter_seq);
    return SUCCESS;
}
