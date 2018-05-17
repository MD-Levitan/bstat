#ifndef __BSTAT_GENERATION_H
#define __BSTAT_GENERATION_H

#include <bstat.h>
#include <cmm.h>
#include "cmms.h"
#include "hmm.h"


// Generation from lib of model TO SEQUENCE
byte generate_cmm_sequence(sequence *seq, cmm_model *model);
byte generate_cmms_sequence(sequence *seq, cmms_model *model);
byte generate_hmm_sequence(sequence *seq, hmm_model *model);
byte generate_dcmm_sequence(sequence *seq, dcmm_model *model);


// Generation TO FILE(ostream)
byte generate_cmm_ostream(ostream *ctx, cmm_model *model, qword T);
byte generate_cmms_ostream(ostream *ctx, cmms_model *model, qword T);
byte generate_hmm_ostream(ostream *ctx, hmm_model *model, qword T);
byte generate_dcmm_ostream(ostream *ctx, dcmm_model *model, qword T);

//
byte generate_to_file(char *filename, void* model, byte model_type, byte isASCII, qword T);

#endif //__BSTAT_GENERATION_H
