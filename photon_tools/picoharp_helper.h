#pragma once

extern "C" {

struct timestamps picoharp_read_file(const char* filename,
                                     unsigned int channel);

struct timestamps {
        double jiffy;
        const uint64_t* timestamps;
        unsigned int nrec;
};

} /* extern "C" */
