#include <cstdint>
#include <iostream>
#include <vector>

const double resolution = 1e-7;

uint64_t acorr(std::vector<uint64_t>& v, uint64_t dt) {
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
	const unsigned int chunk_sz = 1024;
	std::vector<uint64_t> timestamps;
	unsigned int i=0;
	do {
		timestamps.resize(i + chunk_sz);
		std::cin.read((char*) &timestamps[i], chunk_sz);
		i += chunk_sz;
	} while (!std::cin.eof() && !std::cin.fail());

	for (uint64_t dt=0; dt < 100; dt++) {
		double value = acorr(timestamps, dt);
		std::cout.write((char*) &dt, sizeof(uint64_t));
		std::cout.write((char*) &value, sizeof(double));
	}

	return 0;
}

