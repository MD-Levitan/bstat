#include <stdio.h>
#include <time.h>
#include "cmms.h"

//////////////////////////////////////////////////////////////////////
///////////       TESTS OF CMMS AND STATISTIC OF CMMS      ///////////
//////////////////////////////////////////////////////////////////////


int main() {
    srand(time(NULL));
    _entropy = rand;

    printf("//////////////////////////////////////////////////////////////////////\n");
    printf("///////////       TESTS OF CMMS AND STATISTIC OF CMMS      ///////////\n");
    printf("//////////////////////////////////////////////////////////////////////\n");

    sequence seq;
    cmms_model model, model1;
    init_cmms_model(&model, 2, 2);
    init_cmms_model(&model1, 2, 2);
    generate_random_cmms_model(&model);
    generate_random_cmms_model(&model1);

    init_sequence(&seq, 1000);
    //generate_random_hmm_seq(&seq, 2);
    generate_cmms_sequence(&seq, &model);

    printf("Sequence: ");
    for (int j = 0; j < seq.T; ++j) {
        printf("%d ", seq.array[j]);
    }
    printf("\n");

    printf("Estimated model:");
    printf("\n");
    printf("Pi:  ");
    for (int j = 0; j < model.N; ++j) {
        printf("%lf ", model.Pi[j]);
    }
    printf("\n");

    printf("P:  ");
    for (byte i = 0; i < pow(model.N, model.S); ++i) {
        for (int j = 0; j < model.N; ++j) {
            printf("%lf ", model.P[i][j]);
        }
        printf("\n\t");
    }
    printf("\n");


//    double est = 0;
//    printf("estimation = %f\n", est);

    MLE_algorithm_s(&seq, &model1, NULL);
    printf("========MLE========\n");
    printf("\n");
    printf("Pi:  ");
    for (int j = 0; j < model1.N; ++j) {
        printf("%lf ", model1.Pi[j]);
    }
    printf("\n");

    printf("P:  ");
    for (byte i = 0; i < pow(model1.N, model1.S); ++i) {
        for (int j = 0; j < model1.N; ++j) {
            printf("%lf ", model1.P[i][j]);
        }
        printf("\n\t");
    }
    printf("\n");

    generate_random_cmms_model(&model1);
    bootstrap_s(&seq, &model1, 100, NULL);
    printf("========BOOTSTRAP========\n");
    printf("\n");
    printf("Pi:  ");
    for (int j = 0; j < model1.N; ++j) {
        printf("%lf ", model1.Pi[j]);
    }
    printf("\n");

    printf("P:  ");
    for (byte i = 0; i < pow(model1.N, model1.S); ++i) {
        for (int j = 0; j < model1.N; ++j) {
            printf("%lf ", model1.P[i][j]);
        }
        printf("\n\t");
    }
    printf("\n");

    generate_random_cmms_model(&model1);
    smoothed_estimators_s(&seq, &model1, 100, 0.005, NULL);
    printf("========SMOOTHED========\n");
    printf("\n");
    printf("Pi:  ");
    for (int j = 0; j < model1.N; ++j) {
        printf("%lf ", model1.Pi[j]);
    }
    printf("\n");

    printf("P:  ");
    for (byte i = 0; i < pow(model1.N, model1.S); ++i) {
        for (int j = 0; j < model1.N; ++j) {
            printf("%lf ", model1.P[i][j]);
        }
        printf("\n\t");
    }
    printf("\n");


//    int start_time = clock();
//    double average_std = 0;
//    int N = 10;
//    for (int k = 0; k < N; ++k) {
//        estimation_model(&seq, &model1, 0.005, NULL);
//        double std_deviation = 0;
//        std_deviation += standart_deviation_matrix(model.P, model1.P, model.N, model.N);
//        std_deviation += standart_deviation_matrix(model.C, model1.C, model.N, model.M);
//        average_std += std_deviation;
//    }
//    average_std /= N;
//    printf("Time: %f \n", (clock() - (double)start_time) / CLOCKS_PER_SEC);
//
//    printf("Std.Dev.: %f \n", average_std);






//    printf("estimation = %f", est);
//
//    printf("\n=========\n");
//    double stat = 0;
//    statistic_likelihood(&seq, &model1, NULL, 0.005, &stat);
//    printf("Likelihood(std) = %f\n", stat);
//    statistic_likelihood(&seq, &model2, NULL, 0.005, &stat);
//    printf("Likelihood(global) = %f\n", stat);

    free_sequence(&seq);
    free_cmms_model(&model);
    free_cmms_model(&model1);

    return 0;
}
