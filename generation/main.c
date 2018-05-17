#include <generation.h>

int main(){
    printf(" G E N E R A T O R    M O D U L E \n");
    printf("Input file: ");
    char filepath[128];
    scanf("%s", filepath);
    printf("\n");

    printf(" G E N E R A T O R    S E L E C T I O N \n");
    printf("[0] CMM model \n");
    printf("[1] CMM-s model \n");
    printf("[0] HMM model \n");
    printf("[0] DCMM model \n");
    byte model_type = 0;
    printf("Enter choice: ");
    scanf("%d", &model_type);
    printf("\n");

    void *model_ptr;
    if(model_type == 0 || model_type > 3){
        cmm_model model;
        printf("Enter parameters of model: \n");
        printf("Input size of the set of states: ");
        byte n;
        scanf("%d", &n);
        printf("\n");
        init_cmm_model(&model, n);
        printf("Input vector Pi: ");
        printf("\n");
        for(byte i = 0; i < n; ++i){
            scanf("%f", &model.Pi[i]);
        }

        printf("Input matrix P: ");
        for(byte i = 0; i < n; ++i){
            printf("\n");
            for(byte j = 0; j < n; ++j) {
                scanf("%f", &model.P[i][j]);
            }
        }
        model_ptr = &model;
    }

    if(model_type == 1){
        cmms_model model;
        printf("Enter parameters of model: \n");
        printf("Input size of the set of states: ");
        byte n;
        scanf("%d", &n);
        printf("\n");
        printf("Input order: ");
        byte order;
        scanf("%d", &order);
        printf("\n");
        init_cmms_model(&model, n, order);
        printf("Input vector Pi: ");
        printf("\n");
        for(byte i = 0; i < n; ++i){
            scanf("%f", &model.Pi[i]);
        }

        printf("Input matrix P: ");
        for(byte i = 0; i < pow(n, order); ++i){
            printf("\n");
            for(byte j = 0; j < n; ++j) {
                scanf("%f", &model.P[i][j]);
            }
        }
        model_ptr = &model;
    }

    if(model_type == 2){
        hmm_model model;
        printf("Enter parameters of model: \n");
        printf("Input size of the hidden set of states: ");
        byte n;
        scanf("%d", &n);
        printf("\n");
        printf("Input size of the visible set of states: ");
        byte m;
        scanf("%d", &m);
        printf("\n");
        init_hmm_model(&model, n, m);
        printf("Input vector Pi: ");
        printf("\n");
        for(byte i = 0; i < n; ++i){
            scanf("%f", &model.Pi[i]);
        }

        printf("Input matrix P: ");
        for(byte i = 0; i < n; ++i){
            printf("\n");
            for(byte j = 0; j < n; ++j) {
                scanf("%f", &model.P[i][j]);
            }
        }

        printf("Input matrix C: ");
        for(byte i = 0; i < n; ++i){
            printf("\n");
            for(byte j = 0; j < m; ++j) {
                scanf("%f", &model.C[i][j]);
            }
        }

        model_ptr = &model;
    }

    if(model_type == 3){
        dcmm_model model;
        printf("Enter parameters of model: \n");
        printf("Input size of the hidden set of states: ");
        byte n;
        scanf("%d", &n);
        printf("\n");
        printf("Input size of the visible set of states: ");
        byte m;
        scanf("%d", &m);
        printf("\n");
        init_dcmm_model(&model, n, m);
        printf("Input vector Pi: ");
        printf("\n");
        for(byte i = 0; i < n; ++i){
            scanf("%f", &model.Pi[i]);
        }

        printf("Input matrix P: ");
        for(byte i = 0; i < n; ++i){
            printf("\n");
            for(byte j = 0; j < n; ++j) {
                scanf("%f", &model.P[i][j]);
            }
        }

        printf("Input matrix C: ");
        for(byte k = 0; k < n; ++k){
            for(byte i = 0; i < m; ++i) {
                printf("\n");
                for (byte j = 0; j < m; ++j) {
                    scanf("%f", &model.C[k][i][j]);
                }
            }
        }

        model_ptr = &model;
    }

    printf("Input File Format:\n");
    printf("[0] ASCII - A sequence of ASCII 0's and 1's\n");
    printf("[1] Binary - Each byte in data file contains 8 bits of data\n");
    byte isASCII = 0;
    printf("Enter choice:");
    scanf("%d", &isASCII);
    printf("\n");

    printf("Input size of sequence: ");
    qword T = 0;
    scanf("%lli", &T);
    printf("\n");
    byte rv = generate_to_file(filepath, model_ptr, model_type, isASCII, T);
    switch (model_type) {
        case 0:
            free_cmm_model(model_ptr);
            break;

        case 1:
            free_cmms_model(model_ptr);
            break;

        case 2:
            free_hmm_model(model_ptr);
            break;

        case 3:
            free_dcmm_model(model_ptr);
            break;

        default:
            free_cmm_model(model_ptr);
            break;
    }
}