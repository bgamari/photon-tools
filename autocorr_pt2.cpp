#include <vector>
#include <fstream>
#include "pt2.h"

const double resolution = 1e-7;

uint64_t acorr(std::vector<uint64_t> v, uint64_t dt) {
	std::vector<uint64_t>::iterator iter1 = v.begin();
	std::vector<uint64_t>::iterator iter2 = v.begin() + dt;
	uint64_t count;

	while (iter1 != v.end() && iter2 != v.end()) {
		uint64_t diff = *iter2 - *iter1;
		if (diff == 0) count++;
		if (diff <= 0) iter1++;
		if (diff >= 0) iter2++;
	}
	return count;
}

int main(int argc, char** argv) {
	std::vector<std::vector<uint64_t> > chans;
	std::ifstream is(argv[1]);
	pt2_file pt2(is);
	double scale = PT2_TIME_UNIT / resolution;
	unsigned int n_rec = pt2.tttr_hdr.n_records;

	for (int i=0; i < n_rec; i++) {
		pt2_record rec = pt2.read_record();
		if (!rec.special) {
			uint64_t time = (uint64_t) rec.time * scale;
			chans[rec.channel].push_back(time);
		}
	}

	for (uint64_t dt=0; dt < 100; dt++)
		printf("%llu	%llu\n", dt, acorr(ch1, dt));

	return 0;
}

