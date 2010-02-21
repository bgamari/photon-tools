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

struct daton {
  uint64_t abscissa;
  int ordinate;  
};

/*
 * acorr: Compute autocorrelation function a lag given by dt
 *
 * Data abscissae must be monotonic.
 *
 */
uint64_t acorr(std::vector<daton>& v, uint64_t dt) {
	std::vector<daton>::iterator iter1 = v.begin(), iter2 = v.begin();
	uint64_t res = 0;

        // Skip to beginning of shifted vector
        while (iter2->abscissa < dt)
                iter2++;

	while (iter1 != v.end() && iter2 != v.end()) {
                uint64_t diff = iter1->abscissa - (iter2->abscissa - dt);
		if (diff == 0) res += iter1->ordinate * iter2->ordinate;
		if (diff <= 0) iter1++;
		if (diff >= 0) iter2++;
	}
	return res;
}

int main(int argc, char** argv) {
	std::vector<daton> data;
        int min=0, step=1, max=100;

	if (argc == 4) {
		min = atoi(argv[1]);
		max = atoi(argv[2]);
		step = atoi(argv[3]);
	} else {
		printf("Usage: %s MIN MAX STEP\n", argv[0]);
		return 1;
	}

	do {
                daton d;
		std::cin.read((char*) &d.abscissa, sizeof(uint64_t));
                std::cin.read((char*) &d.ordinate, sizeof(int));
                data.push_back(d);
	} while (!std::cin.eof() && !std::cin.fail());

	for (uint64_t dt=min; dt < max; dt += step) {
		double value = acorr(data, dt);
		std::cout.write((char*) &dt, sizeof(uint64_t));
		std::cout.write((char*) &value, sizeof(double));
                std::cerr << dt << "\t" << value << "\n";
	}

	return 0;
}

