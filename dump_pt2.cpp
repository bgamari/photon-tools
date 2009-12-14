#include <fstream>
#include <cstdio>
#include "pt2.h"

int main(int argc, char** argv) {
	std::ifstream is(argv[1]);
	pt2_file pt2(is);

	for (int i=0; i < pt2.tttr_hdr.n_records; i++) {
		pt2_record rec = pt2.read_record();
		if (rec.special) {
			printf("# ");
			if (rec.markers == 0)
				printf("overflow\n");
			else
				printf("marker %d\n", rec.markers);
		} else
			printf("%7llu\t%d\n", rec.time, rec.channel);
	}
	return 0;
}

