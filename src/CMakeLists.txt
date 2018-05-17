set(src
        cmm.c
        cmms.c
        dcmm.c
        hmm.c
        statistic.c
        util.c
        tests.c
        sequence.c processing.c stream.c)


add_library(bstat_src OBJECT ${src})
set_property(TARGET bstat_src PROPERTY POSITION_INDEPENDENT_CODE)
set_property(TARGET bstat_src PROPERTY POSITION_INDEPENDENT_CODE ON)