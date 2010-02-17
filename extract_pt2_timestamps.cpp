#include <iostream>
#include "pt2.h"

const double resolution = 1e-7;

int main(int argc, char** argv) {
	pt2_file pt2(std::cin);
	double scale = PT2_TIME_UNIT / resolution;
	unsigned int n_rec = pt2.tttr_hdr.n_records;

	for (int i=0; i < n_rec; i++) {
		pt2_record rec = pt2.read_record();
		if (!rec.special) {
			uint64_t time = (uint64_t) rec.time * scale;
			std::cout.write((char*) &time, sizeof(uint64_t));
		}
	}

	return 0;
}

