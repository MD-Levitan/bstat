#include "tests.h"


double *statistic_cmms(byte N, byte S, qword len, byte alg){
    
    
    sequence seq;
    init_sequence(&seq, len);
    generate_random_sequence(&seq, N);
    cmms_model model;
    init_cmms_model(&model, N, S);
    double *n;
    _memcheck(n, pow(N, S) * sizeof(double));
    switch(alg){
        case 0:
            MLE_algorithm_s(&seq, &model, n);
            break;
        case 1:
            bootstrap_s(&seq, &model, 1000, n);
            break;
        case 2:
            smoothed_estimators_s(&seq, &model, 1000, 0.5, n);
            break;
        default:
            MLE_algorithm_s(&seq, &model, n);
            break;
    }
    double stat;
    statistic_chi_square_cmms(N, S, n, &model, NULL, 0.005, &stat);
    byte _m = 0;
    
    for(byte i = 0; i < pow(model.N, model.S); ++i){
        for(byte j = 0; j < model.N; ++j) {
            if (model.P[i][j] != 0)
                _m++;
        }
    }
    double *res;
    _memcheck(res, 2 * sizeof(double));
    res[0] = stat;
    res[1] = _m;
    free_sequence(&seq);
    free_cmms_model(&model);
    free(n);
    return res;
}

double *statistic_cmms_2(byte N, byte S, qword len, byte alg){
    cmms_model model_bad;
    init_cmms_model(&model_bad, N, S);
    generate_uniform_cmms_model(&model_bad);
    model_bad.P[0][0] = 0.4;
    model_bad.P[0][1] = 0.6;
    model_bad.P[1][0] = 0.65;
    model_bad.P[1][1] = 0.35;
    model_bad.P[2][0] = 0.6;
    model_bad.P[2][1] = 0.4;
    model_bad.P[3][0] = 0.35;
    model_bad.P[3][1] = 0.65;
    
    
    
    sequence seq;
    init_sequence(&seq, len);
    generate_cmms_sequence(&seq, &model_bad);
    free_cmms_model(&model_bad);
    
    cmms_model model;
    init_cmms_model(&model, N, S);
    double *n;
    _memcheck(n, pow(N, S) * sizeof(double));
    switch(alg){
        case 0:
            MLE_algorithm_s(&seq, &model, n);
            break;
        case 1:
            bootstrap_s(&seq, &model, 1000, n);
            break;
        case 2:
            smoothed_estimators_s(&seq, &model, 1000, 0.5, n);
            break;
        default:
            MLE_algorithm_s(&seq, &model, n);
            break;
    }
    double stat;
    statistic_chi_square_cmms(N, S, n, &model, NULL, 0.005, &stat);
    byte _m = 0;
    for(byte i = 0; i < pow(model.N, model.S); ++i){
        for(byte j = 0; j < model.N; ++j) {
            if (model.P[i][j] != 0)
                _m++;
        }
    }
    double *res;
    _memcheck(res, 2 * sizeof(double));
    res[0] = stat;
    res[1] = _m;
    free_sequence(&seq);
    free_cmms_model(&model);
    
    free(n);
    return res;
}

double statistic_cmms_3(byte N, byte S, qword len, byte alg){
    
    time_t start = clock();
    
    sequence seq;
    init_sequence(&seq, len);
    generate_random_sequence(&seq, N);
    cmms_model model;
    init_cmms_model(&model, N, S);
    double *n;
    _memcheck(n, pow(N, S) * sizeof(double));
    switch(alg){
        case 0:
            MLE_algorithm_s(&seq, &model, n);
            break;
        case 1:
            bootstrap_s(&seq, &model, 1000, n);
            break;
        case 2:
            smoothed_estimators_s(&seq, &model, 1000, 0.5, n);
            break;
        default:
            MLE_algorithm_s(&seq, &model, n);
            break;
    }
    double stat;
    statistic_chi_square_cmms(N, S, n, &model, NULL, 0.005, &stat);
    byte _m = 0;
    for(byte i = 0; i < pow(model.N, model.S); ++i){
        for(byte j = 0; j < model.N; ++j) {
            if (model.P[i][j] != 0)
                _m++;
        }
    }
    free_sequence(&seq);
    free_cmms_model(&model);
    free(n);
    return  (double)(clock() - start) / CLOCKS_PER_SEC;
}

double statistic_dcmm(byte N, byte M, qword len){
    sequence seq;
    init_sequence(&seq, len);
    generate_random_sequence(&seq, M);
    dcmm_model model;
    init_dcmm_model(&model, N, M);
    estimation_model_d(&seq, &model, 0.005, NULL);
    double stat;
   // statistic_likelihood_dcmm(&seq, &model, NULL, 0.005, &stat);
    free_sequence(&seq);
    free_dcmm_model(&model);
    return  stat;
}

double statistic_dcmm_3(byte N, byte M, qword len){
    
    time_t start = clock();
    
    sequence seq;
    init_sequence(&seq, len);
    generate_random_sequence(&seq, M);
    dcmm_model model;
    init_dcmm_model(&model, N, M);
    estimation_model_d(&seq, &model, 0.005, NULL);
    double stat;
    //statistic_likelihood_dcmm(&seq, &model, NULL, 0.005, &stat);
    free_sequence(&seq);
    free_dcmm_model(&model);
    return  (clock() - (double)start) / CLOCKS_PER_SEC;
}


double *statistic_cmm(byte N, qword len, byte alg){
    sequence seq;
    read_file_2_sequence(&seq, "sampledata-100MB.bin", N);
    //init_sequence(&seq, len);
    //generate_random_sequence(&seq, N);
    
    cmm_model model;
    init_cmm_model(&model, N);
    double *n;
    _memcheck(n, N * sizeof(double));
    switch(alg){
        case 0:
            MLE_algorithm(&seq, &model, n);
            break;
        case 1:
            bootstrap(&seq, &model, 1000, n);
            break;
        case 2:
            smoothed_estimators(&seq, &model, 1000, 0.5, n);
            break;
        default:
            MLE_algorithm(&seq, &model, n);
            break;
    }
    double stat;
    statistic_chi_square_cmm(N, n, &model, NULL, 0.005, &stat);
    byte _m = 0;
    for(byte i = 0; i < model.N; ++i){
        for(byte j = 0; j < model.N; ++j) {
            if (model.P[i][j] != 0)
                _m ++;
        }
    }
    double *res;
    _memcheck(res, 2 * sizeof(double));
    res[0] = stat;
    res[1] = _m;
    free_sequence(&seq);
    free_cmm_model(&model);
    free(n);
    return res;
}
