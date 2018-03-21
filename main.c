#include <stdio.h>
#include <time.h>
#include "hmm.h"
#include "statistic.h"

int main() {
    srand(time(NULL));
    _entropy = rand;
    hmm_seq seq;
    hmm_model model, model1, model2;
    init_hmm_model(&model, 2, 2);
    init_hmm_model(&model1, 2, 2);
    init_hmm_model(&model2, 2, 2);
    generate_random_hmm_model(&model);
    generate_random_hmm_model(&model1);
    generate_random_hmm_model(&model2);

    init_hmm_seq(&seq, 1000, &model);
    //generate_random_hmm_seq(&seq, 2);
    generate_hmm_seq(&seq, &model);
    double est = 0;
    estimation_model(&seq, &model1, 0.005, &est);
    printf("estimation = %f\n", est);
    estimation_model_gl(&seq, &model2, 4, 0.005, &est);

//    printf("\n Seq: ");
//    for (int j = 0; j < seq.T; ++j) {
//        printf("%d ", seq.array[j]);
//    }
//    printf("\n");

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


    printf("Pi:  ");
    for (int j = 0; j < model.N; ++j) {
        printf("%lf ", model.Pi[j]);
    }

    printf("\nP:  ");
    for (byte i = 0; i < model.N; ++i) {
        for (int j = 0; j < model.N; ++j) {
            printf("%lf ", model.P[i][j]);
        }
        printf("\n\t");
    }

    printf("\nC:  ");
    for (byte i = 0; i < model.N; ++i) {
        for (int j = 0; j < model.M; ++j) {
            printf("%lf ", model.C[i][j]);
        }
        printf("\n\t");
    }




    printf("NEW MODEL\n");
    printf("Pi:  ");
    for (int j = 0; j < model1.N; ++j) {
        printf("%lf ", model1.Pi[j]);
    }

    printf("\nP:  ");
    for (byte i = 0; i < model1.N; ++i) {
        for (int j = 0; j < model1.N; ++j) {
            printf("%lf ", model1.P[i][j]);
        }
        printf("\n\t");
    }

    printf("\nC:  ");
    for (byte i = 0; i < model1.N; ++i) {
        for (int j = 0; j < model1.M; ++j) {
            printf("%lf ", model1.C[i][j]);
        }
        printf("\n\t");
    }
    printf("estimation = %f", est);
    
    printf("\n=========\n");
    double stat = 0;
    statistic_likelihood(&seq, &model1, NULL, 0.005, &stat);
    printf("Likelihood(std) = %f\n", stat);
    statistic_likelihood(&seq, &model2, NULL, 0.005, &stat);
    printf("Likelihood(global) = %f\n", stat);

    free_hmm_seq(&seq);
    free_hmm_model(&model);
    free_hmm_model(&model1);

    return 0;
}