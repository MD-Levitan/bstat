#include "bstat.h"
byte f_open(stream *ctx, const char* filename, byte m){
    if (!ctx) {
        return ERROR;
    }
    if (m % 2){
        perror("Error in size param");
        return ERROR;
    }
    ctx->m = m;
    ctx->bits = (byte)log2(ctx->m);
    ctx->file = fopen(filename, "r");
    if(ctx->file) {
        perror("File opening failed");
        return ERROR;
    }
}

void f_close(stream *ctx){
    if(!ctx || !ctx->file)
        return;
    fclose(ctx->file);
}

byte is_open(stream *ctx){
    if(ctx && ctx->file)
        return SUCCESS;
    return ERROR;
}

byte is_end(stream *ctx){
    if(ctx->buffer_size < ctx->bits) {
        byte rv = 0;
        byte _buffer = ctx->buffer;
        rv = (byte)fread(&ctx->buffer, 1, 1, ctx->file);
        ctx->buffer <<= ctx->buffer_size;
        ctx->buffer += _buffer;
        if(!rv)
            return 1;
        ctx->buffer_size += 8;
    }
    return 0;
}

byte get(stream *ctx) {
    if (ctx->buffer_size < ctx->bits) {
        byte rv = 0;
        byte _buffer = ctx->buffer;
        rv = (byte)fread(&ctx->buffer, 1, 1, ctx->file);
        ctx->buffer <<= ctx->buffer_size;
        ctx->buffer += _buffer;
        ctx->buffer_size += 8;
    }
    byte res = ctx->buffer % ctx->m;
    ctx->buffer >>= ctx->bits;
    ctx->buffer_size -= ctx->bits;
    return res;
}

qword get_size(stream *ctx){
    if(!is_open(ctx))
        return 0;
    fseek(ctx->file, 0L, SEEK_END);
    qword size = ftell(ctx->file);
    fseek(ctx->file, 0L, SEEK_SET);
    return size * (1 << 8 / ctx->m);
}