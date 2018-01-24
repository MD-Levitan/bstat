#include <stdio.h>
#include <time.h>
#include "hmm.h"

int main() {
    srand(time(NULL));
    _entropy = rand;
    hmm_seq seq;
    hmm_model model, model1;
    init_hmm_model(&model, 2, 2);
    init_hmm_model(&model1, 2, 2);
    generate_random_hmm_model(&model);
    generate_random_hmm_model(&model1);

    init_hmm_seq(&seq, 1000, &model);
    //generate_random_hmm_seq(&seq, 2);
    generate_hmm_seq(&seq, &model);

    int start_time = clock();
    estimation_model(&seq, &model1, 0.000000001);
    printf("Time: %f \n", (clock() - (double)start_time) / CLOCKS_PER_SEC);

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

    printf("\n Seq: ");
    for (int j = 0; j < seq.T; ++j) {
        printf("%d ", seq.array[j]);
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


    free_hmm_seq(&seq);
    free_hmm_model(&model);
    free_hmm_model(&model1);

    return 0;
}