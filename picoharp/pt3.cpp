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

#include "pt3.h"
#include <cstring>
#include <cstdlib>
#include <fstream>

pt3_record pt3_file::read_record() {
	uint32_t rec;
	pt3_record r;

	records_read++;
	if (records_read > tttr_hdr.n_records)
		throw std::runtime_error("Read too many records");

	is.read((char*) &rec, 4);
        r.numsync = 0xffff & rec;
	r.channel = 0xf & (rec >> 28);
	
	if (r.channel == 0xf) {
		r.is_special = true;
                r.special.markers = 0x0fff & (rec >> 16);
		if (r.special.markers == 0) { // overflow record
			overflow_time += PT3_WRAPAROUND_TIME;
                }
	} else {
                r.is_special = false;
                r.normal.time = 0x0fff & (rec >> 16);
                r.normal.time += (overflow_time + r.numsync) * sync_period;
        }

	return r;
}

uint64_t* pt3_file::get_timestamps(unsigned int channel,
                                   unsigned int *n_records)
{
    unsigned int n_rec = tttr_hdr.n_records;
    uint64_t *buffer = new uint64_t[n_rec];
    unsigned int n = 0;

    for (int i=0; i < n_rec; i++) {
	pt3_record rec = read_record();
	if (!rec.is_special && (channel == 0xf || channel == rec.channel)) {
	    buffer[i] = rec.normal.time;
	    n++;
	}
    }
    
    uint64_t *new_buffer = (uint64_t *) malloc(n * sizeof(uint64_t));
    memcpy(new_buffer, buffer, n*sizeof(uint64_t));
    delete[] buffer;
    *n_records = n;
    return new_buffer;
}

