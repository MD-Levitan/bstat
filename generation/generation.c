#include <generation.h>

byte generate_cmm_ostream(ostream *ctx, cmm_model *model, qword T){
    if(_entropy == NULL || iopened(ctx) || model == NULL || T < 1)
        return ERROR;
    byte rv = 0;
    rv = osetm(ctx, model->N);
    if(rv != SUCCESS)
        return rv;
    byte visible = ch(model->Pi, model->N);
    for (qword i = 0; i < T; ++i) {
        oput(ctx, visible);
        visible = ch(model->P[visible], model->N);
    }
}

byte generate_cmms_ostream(ostream *ctx, cmms_model *model, qword T){
    if(_entropy == NULL || iopened(ctx) || model == NULL || T < (model->S / 8))
        return ERROR;
    byte rv = 0;
    rv = osetm(ctx, model->N);
    if(rv != SUCCESS)
        return rv;

    qword current_word = 0;
    byte visible;
    qword num_bit = pow(model->N, model->S - 1);
    for (qword i = 0; i < model->S; ++i) { //*check
        visible = ch(model->Pi, model->N);
        current_word *= model->N;
        current_word += visible;
        oput(ctx, visible);
    }
    for (qword i = model->S; i < T; ++i) {
        visible = ch(model->P[current_word], model->N);
        oput(ctx, visible);
        current_word %= num_bit;
        current_word *= model->N;
        current_word += visible;
    }
}

byte generate_hmm_ostream(ostream *ctx, hmm_model *model, qword T){
    if(_entropy == NULL || iopened(ctx) || model == NULL || T < 1)
        return ERROR;
    byte rv = 0;
    rv = osetm(ctx, model->M);
    if(rv != SUCCESS)
        return rv;

    byte hidden = ch(model->Pi, model->N);
    byte visible = ch(model->C[hidden], model->M);
    for (qword i = 0; i < T; ++i) {
        oput(ctx, visible);
        hidden = ch(model->P[hidden], model->N);
        visible = ch(model->C[hidden], model->M);
    }
}

byte generate_dcmm_ostream(ostream *ctx, dcmm_model *model, qword T){
    if(_entropy == NULL || iopened(ctx) || model == NULL || T < 1)
        return ERROR;
    byte rv = 0;
    rv = osetm(ctx, model->M);
    if(rv != SUCCESS)
        return rv;

    double hidden_start[model->M];
    double var = 1.0 / model->M;
    for(byte i = 0; i < model->M; ++i){
        hidden_start[i] = var;
    }
    byte hidden = ch(model->Pi, model->N);
    byte visible = ch(hidden_start, model->M);
    visible = ch(model->C[hidden][visible], model->M);
    for (qword i = 0; i < T; ++i) {
        oput(ctx, visible);
        hidden = ch(model->P[hidden], model->N);
        visible = ch(model->C[hidden][visible], model->M);
    }
}

byte generate_to_file(char *filename, void* model, byte model_type, byte isASCII, qword T){
    ostream os;
    int rv = 0;
    rv = oopen(&os, filename, isASCII);
    if(rv != SUCCESS){
        printf(" ERROR in oopen file");
    }
    switch (model_type) {
        case 0:
            rv = generate_cmm_ostream(&os, (cmm_model *) model, T);
            break;

        case 1:
            rv = generate_cmms_ostream(&os, (cmms_model *) model, T);
            break;

        case 2:
            rv = generate_hmm_ostream(&os, (hmm_model *) model, T);
            break;

        case 3:
            rv = generate_dcmm_ostream(&os, (dcmm_model *) model, T);
            break;

        default:
            rv = generate_cmm_ostream(&os, (cmm_model *) model, T);
            break;
    }
    oclose(&os);
    return rv;
}