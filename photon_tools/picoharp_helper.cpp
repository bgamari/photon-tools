#include "picoharp.h"
#include "picoharp_helper.h"
#include <istream>
#include <fstream>

extern "C" {

struct timestamps picoharp_read_file(const char* filename,
                                     unsigned int channel)
{
        std::ifstream is(filename);
        picoharp_file pf(is);
        unsigned int nrec = 0;
        const uint64_t* data = pf.get_timestamps(channel, &nrec);
        const double jiffy = 1e-9 * pf.board_hdr.resolution;
        const timestamps ret = {jiffy, data, nrec};
        return ret;
}

} /* extern "C" */
