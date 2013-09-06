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


#include "pt2.h"
#include <cstring>
#include <cstdlib>
#include <fstream>

pt2_record pt2_file::read_record() {
	uint32_t rec;
	pt2_record r;

	records_read++;
	if (records_read > tttr_hdr.n_records)
		throw std::runtime_error("Read too many records");

	is.read((char*) &rec, 4);
	r.time = 0x0fffffff & rec;
	r.channel = (0xf0000000 & rec) >> 28;
	
	if (r.channel == 0xf) {
		r.special = true;
                r.time &= ~0xf; // the bottom 4 bits are markers in this case
		r.markers = r.time & 0xfLL;
		if (r.markers == 0) // overflow record
			overflow_time += PT2_WRAPAROUND_TIME;
	} else {
                r.special = false;
                r.markers = 0;
        }
	r.time += overflow_time;

	return r;
}

std::vector<pt2_record> pt2_file::read_all_records() {
	int n = tttr_hdr.n_records;
	std::vector<pt2_record> records(n);
	for (int i=0; i < n; i++)
		records.push_back(read_record());
	return records;
}

uint64_t* pt2_file::get_timestamps(unsigned int channel,
                                   unsigned int *n_records)
{
    unsigned int n_rec = tttr_hdr.n_records;
    uint64_t *buffer = new uint64_t[n_rec];
    unsigned int n = 0;

    for (int i=0; i < n_rec; i++) {
	pt2_record rec = read_record();
	if (!rec.special && (channel == 0xf || channel == rec.channel)) {
	    buffer[i] = rec.time;
	    n++;
	}
    }
    
    uint64_t *new_buffer = (uint64_t *) malloc(n * sizeof(uint64_t));
    memcpy(new_buffer, buffer, n*sizeof(uint64_t));
    delete[] buffer;
    *n_records = n;
    return new_buffer;
}

