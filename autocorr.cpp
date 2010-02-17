/* fcs-tools - Tools for FCS data analysis
 *
 * Copyright © 2010 Ben Gamari
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/ .
 *
 * Author: Ben Gamari <bgamari@physics.umass.edu>
 */


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

