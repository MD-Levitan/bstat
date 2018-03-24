#include <libnet.h>
#include <hmm.h>
#include "cmms.h"

/////////////////////////////////////////////////
//////////////// INIT, FREE, COPY ///////////////
/////////////////////////////////////////////////

void init_cmms_model(cmms_model *ctx, byte n, byte s){
    if(!ctx)
        return;
    ctx->N = n;
    ctx->S = s;
    word ns = pow(ctx->N, s);
    _memcheck(ctx->P, ns * sizeof(double*));
    _memcheck(ctx->Pi, ctx->N * sizeof(double));

    for (byte i = 0; i < ns; ++i)
        _memcheck(ctx->P[i], ctx->N * sizeof(double));
}

void init_cmms_sequence(sequence *ctx, qword t){
    if(!ctx)
        return;

    ctx->T = t;
    _memcheck(ctx->array, ctx->T * sizeof(byte));
}

void copy_cmms_model(cmms_model *dist, cmms_model *src){
    if(!dist || !src || dist->N != src->N || dist->S != src->S)
        return;

    for (int j = 0; j < src->N; ++j) {
        dist->Pi[j] = src->Pi[j];
    }

    for (byte i = 0; i < pow(src->N, src->S); ++i) {
        for (int j = 0; j < src->N; ++j) {
            dist->P[i][j] = src->P[i][j];
        }
    }
}

void free_cmms_model(cmms_model *ctx){
    if(!ctx)
        return;

    for (byte i = 0; i < pow(ctx->N, ctx->S); ++i)
        free(ctx->P[i]);

    free(ctx->Pi);
    free(ctx->P);

    ctx->N = 0;
    ctx->S = 0;
}

/////////////////////////////////////////////////


byte generate_cmms_sequence(sequence *seq, cmms_model *model){
    if(_entropy == NULL || seq == NULL || model == NULL || seq->T <= model->S)
        return ERROR;

    seq->m = model->N;

    qword current_word = 0;
    byte visible;
    qword num_bit = ((qword) -1) ^ (1 << (model->S - 1));
    for (qword i = 0; i < model->S; ++i) { //*check
        visible = ch(model->Pi, model->N);
        current_word <<= 1;
        current_word += visible;
        seq->array[i] = visible;
    }
    for (qword i = model->S; i < seq->T; ++i) {
        visible = ch(model->P[current_word], model->N);
        seq->array[i] = visible;
        current_word &= num_bit;
        current_word <<= 1;
        current_word += visible;
    }
}

byte generate_cmms_model(cmms_model *model, byte type){
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
        for (byte i = 0; i < pow(model->N, model->S); ++i) {
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
        for (byte i = 0; i < pow(model->N, model->S); ++i) {
            for (byte j = 0; j < model->N; ++j)
                model->P[i][j] = (double) 1 / model->N;
        }
    }

    return SUCCESS;
}

byte MLE_algorithm_s(sequence *seq, cmms_model *model, double *_n) {
    if (!model || !seq || (seq->m && seq->m != model->N) || seq->T <= model->S)
        return ERROR;

    //calculate Pi
    double value = 1.0 / model->N;
    for (byte i = 0; i < model->N; ++i) {
        model->Pi[i] = value;
    }

    qword num_bit = ((qword) -1) ^ (1 << (model->S - 1));
    qword ns = pow(model->N, model->S);
    qword *n;
    _memcheck(n, ns * sizeof(qword));
    for (byte i = 0; i < ns; ++i) {
        n[i] = 0;
    }

    for (byte i = 0; i < ns; ++i) {
        for (byte j = 0; j < model->N; ++j) {
            model->P[i][j] = 0;
        }
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

    //calculate n_i and P
    //Complexity = O(2*N*N) + 0(T-1)


    qword current_word = 0;
    for (qword t = 0; t < model->S - 1; ++t) {
        if(seq->array[t] < model->N){
            current_word <<= 1;
            current_word += seq->array[t];
        }else{
            return ERROR;
        }
    }
    for (qword t = model->S - 1; t < seq->T - 1; ++t) {
        current_word <<= 1;
        current_word += seq->array[t];
        if(seq->array[t] < model->N){
            n[current_word]++;
            model->P[current_word][seq->array[t + 1]]++;
        }else{
            return ERROR;
        }
        current_word &= num_bit;
    }


    for (byte i = 0; i < ns; ++i) {
        for (byte j = 0; j < model->N; ++j) {
            model->P[i][j] = n[i] ? model->P[i][j] / n[i] : 1.0 / model->N;
        }
    }
    if(_n != NULL) {
        for (byte i = 0; i < ns; ++i) {
            _n[i] = (double) n[i];
        }
    }
    free(n);
    return SUCCESS;
}

byte bootstrap_s(sequence *seq, cmms_model *model, word repeats, double *_n){
    if (!model || !seq || (seq->m && seq->m != model->N))
        return ERROR;
    cmms_model iter_model;
    sequence iter_seq;
    cmms_model average_model;
    init_cmms_model(&iter_model, model->N, model->S);
    init_cmms_model(&average_model, model->N, model->S);
    double *n, *n_average;
    qword ns = pow(model->N, model->S);
    _memcheck(n, ns * sizeof(double));
    _memcheck(n_average, ns * sizeof(double));


    _check(MLE_algorithm_s(seq, &iter_model, n));
    memcpy(n_average, n, model->N * sizeof(double));
    copy_cmms_model(&average_model, &iter_model);
    init_cmms_sequence(&iter_seq, seq->T);
    for (word it = 1; it < repeats; ++it) {
        generate_cmms_sequence(&iter_seq, &iter_model);
        _check(MLE_algorithm_s(&iter_seq, &iter_model, n));
        for(byte i = 0; i < model->N; ++i){
            n_average[i] += n[i];
        }
        for(byte i = 0; i < ns; ++i){
            for(byte j = 0; j < model->N; ++j){
                average_model.P[i][j] += iter_model.P[i][j];
            }
        }
    }
    for(byte i = 0; i < ns; ++i){
        for(byte j = 0; j < model->N; ++j){
            average_model.P[i][j] /= repeats;
        }
    }
    for(byte i = 0; i < model->N; ++i){
        n_average[i] /= repeats;
    }
    if(_n)
        memcpy(_n, n, model->N * sizeof(double));
    free(n);
    free(n_average);
    copy_cmms_model(model, &average_model);
    free_cmms_model(&average_model);
    free_cmms_model(&iter_model);
    free_sequence(&iter_seq);
    return SUCCESS;
}

byte smoothed_estimators_s(sequence *seq, cmms_model *model, word repeats, double u, double *_n) {
    if (!model || !seq || (seq->m && seq->m != model->N))
        return ERROR;
    cmms_model iter_model;
    sequence iter_seq;
    cmms_model average_model;
    init_cmms_model(&iter_model, model->N, model->S);
    init_cmms_model(&average_model, model->N, model->S);
    double gamma = pow(seq->T,-u);
    double omega = 1 + gamma * model->N;
    double *n, *n_average;
    qword ns = pow(model->N, model->S);
    _memcheck(n, ns * sizeof(double));
    _memcheck(n_average, ns * sizeof(double));


    _check(MLE_algorithm_s(seq, &iter_model, n));
    copy_cmms_model(&average_model, &iter_model);
    for(byte i = 0; i < ns; ++i){
        for(byte j = 0; j < model->N; ++j){
            iter_model.P[i][j] = (iter_model.P[i][j] + gamma) / omega;
        }
    }
    init_sequence(&iter_seq, seq->T);
    for (word it = 1; it < repeats; ++it) {
        generate_cmms_sequence(&iter_seq, &iter_model);
        _check(MLE_algorithm_s(&iter_seq, &iter_model, n));
        for(byte i = 0; i < model->N; ++i){
            n_average[i] += n[i];
        }
        for(byte i = 0; i < ns; ++i){
            for(byte j = 0; j < model->N; ++j){
                average_model.P[i][j] += iter_model.P[i][j];
            }
        }
    }
    for(byte i = 0; i < ns; ++i){
        for(byte j = 0; j < model->N; ++j){
            average_model.P[i][j] /= repeats;
        }
    }
    for(byte i = 0; i < model->N; ++i){
        n_average[i] /= repeats;
    }
    if(_n)
        memcpy(_n, n, model->N * sizeof(double));
    free(n);
    free(n_average);
    copy_cmms_model(model, &average_model);
    free_cmms_model(&average_model);
    free_cmms_model(&iter_model);
    free_sequence(&iter_seq);
    return SUCCESS;
}
