/* fcs-tools - Tools for FCS data analysis
 *
 * Copyright © 2013 Ben Gamari
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
#include "pt3.h"

int main(int argc, char** argv) {
	std::ifstream is(argv[1]);
	pt3_file pt3(is);

	printf("# jiffy = %f\n", pt3.board_hdr.resolution);
	for (unsigned int i=0; i < pt3.tttr_hdr.n_records; i++) {
		pt3_record rec = pt3.read_record();
		if (rec.is_special) {
			printf("# ");
			if (rec.special.markers == 0)
				printf("overflow\n");
			else
				printf("marker %d\n", rec.special.markers);
		} else
			printf("%7lu\t%d\n", rec.normal.time, rec.channel);
	}
	return 0;
}

