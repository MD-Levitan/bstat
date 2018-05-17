#include "processing.h"

byte read_file_2_sequence(sequence *seq, const char *filename, byte m) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        return NULL;
    }
    fseek(file, 0L, SEEK_END);
    qword size = ftell(file);
    fseek(file, 0L, SEEK_SET);
    seq->m = m;
    word buffer;
    byte _buffer = 0, buffer_len = 0, len = 0, bits = (byte)log2(seq->m);
    init_sequence(seq, size * 8 / bits);
    qword counter = 0;
    len = fread(&buffer, 1, 1, file);
    while (len){
        buffer <<= buffer_len;
        buffer += _buffer;
        buffer_len += 8;
        while(buffer_len >= bits){
            seq->array[counter++] = buffer % seq->m;
            buffer >>= bits;
            buffer_len -= bits;
        }
        _buffer = buffer;
        len = fread(&buffer, 1, 1, file);
    }
    
    fclose(file);
}
