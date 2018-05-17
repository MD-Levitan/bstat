#include "bstat.h"
byte iopen(istream *ctx, const char* filename, byte m){
    if (!ctx) {
        return ERROR;
    }
    if (m % 2){
        perror("Error in size param");
        return ERROR;
    }
    ctx->buffer_size = 0;
    ctx->m = m;
    ctx->bits = (byte)log2(ctx->m);
    ctx->file = fopen(filename, "r");
    if(ctx->file) {
        perror("File opening failed");
        return ERROR;
    }
    return SUCCESS;

}

void iclose(istream *ctx){
    if(!ctx || !ctx->file)
        return;
    fclose(ctx->file);
}

byte iopened(istream *ctx){
    if(ctx && ctx->file)
        return SUCCESS;
    return ERROR;
}

byte iend(istream *ctx){
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

byte iget(istream *ctx) {
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

qword isize(istream *ctx){
    if(!iopened(ctx))
        return 0;
    fseek(ctx->file, 0L, SEEK_END);
    qword size = ftell(ctx->file);
    fseek(ctx->file, 0L, SEEK_SET);
    return size * (1 << 8 / ctx->m);
}


byte oopen(ostream *ctx, const char* filename, byte isASCII){
    if (!ctx) {
        return ERROR;
    }
    ctx->buffer_size = 0;
    ctx->isASCII = isASCII;
    ctx->file = fopen(filename, "w");
    if(ctx->file) {
        perror("File opening failed");
        return ERROR;
    }
    return SUCCESS;
}

byte osetm(ostream *ctx, byte m){
    if (!ctx) {
        return ERROR;
    }
    if (m % 2){
        perror("Error in size param");
        return ERROR;
    }
    ctx->m = m;
    ctx->bits = (byte)log2(ctx->m);
    return SUCCESS;
}

void  oclose(ostream *ctx){
    if(!ctx || !ctx->file)
        return;
    if(ctx->buffer_size != 0){
        byte buff = ctx->buffer & 0xFF;
        ctx->buffer_size = 0;
        fwrite(&buff, 1, 1, ctx->file);
    }
    fclose(ctx->file);
}

byte  oopened(ostream *ctx){
    if(ctx && ctx->file)
        return SUCCESS;
    return ERROR;
}

void oput_ASCII(ostream *ctx, byte in){
    byte rv = (byte) fwrite(&in, 1, 1, ctx->file);
    return; //rv;
}

void oput_binary(ostream *ctx, byte in){
    ctx->buffer <<= ctx->bits;
    ctx->buffer += in;
    ctx->buffer_size += ctx->m;

    if(ctx->buffer_size >= 8) {
        byte rv;
        byte res = ctx->buffer & 0xFF;
        ctx->buffer >>= 8;
        ctx->buffer_size -= 8;
        rv = (byte) fwrite(&res, 1, 1, ctx->file);
    }
    return; //res;
}

void oput(ostream *ctx, byte in){
    return ctx->isASCII ? oput_ASCII(ctx, in) : oput_binary(ctx, in);
}

