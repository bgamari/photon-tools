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


#include <fstream>
#include <cstdio>
#include "pt2.h"

int main(int argc, char** argv) {
	std::ifstream is(argv[1]);
	pt2_file pt2(is);

	printf("# jiffy = %f\n", pt2.board_hdr.resolution);
	for (unsigned int i=0; i < pt2.tttr_hdr.n_records; i++) {
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

