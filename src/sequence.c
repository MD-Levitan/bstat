#include "bstat.h"

void init_sequence(sequence *ctx, qword t){
    if(!ctx)
        return;
    
    ctx->T = t;
    _memcheck(ctx->array, ctx->T * sizeof(byte));
}

void free_sequence(sequence *ctx){
    if(!ctx)
        return;
    
    free(ctx->array);
    ctx->T = 0;
}