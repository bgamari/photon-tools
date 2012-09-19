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
#include <stdexcept>
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

void pt2_file::read_headers() {
	is.read((char*) &text_hdr, sizeof(pt2_text_hdr));
	if (strcmp(text_hdr.ident, "PicoHarp 300"))
		throw std::runtime_error("File identifier not found");
	if (strncmp(text_hdr.format_version, "2.0", 3))
		throw std::runtime_error("Unsupported file format version");

	is.read((char*) &binary_hdr, sizeof(pt2_binary_hdr));
	if (binary_hdr.meas_mode != PT2_MEASMODE_T2)
		throw std::runtime_error("Unsupported measurement mode");

	is.read((char*) &board_hdr, sizeof(pt2_board_hdr));
	is.read((char*) &tttr_hdr, sizeof(pt2_tttr_hdr));
	is.ignore(4*tttr_hdr.imaging_hdr_sz);
}

uint64_t *get_pt2_timestamps(const char *filename, unsigned int channel, unsigned int *n_records)
{
    std::ifstream is(filename);
    pt2_file *pt2 = new pt2_file(is);
    unsigned int n_rec = pt2->tttr_hdr.n_records;
    uint64_t *buffer = new uint64_t[n_rec];
    unsigned int n = 0;

    for (int i=0; i < n_rec; i++) {
	pt2_record rec = pt2->read_record();
	if (!rec.special && (channel == 0xf || channel == rec.channel)) {
	    buffer[i] = rec.time;
	    n++;
	}
    }
    delete pt2;
    
    uint64_t *new_buffer = (uint64_t *) malloc(n * sizeof(uint64_t));
    memcpy(new_buffer, buffer, n*sizeof(uint64_t));
    delete[] buffer;
    *n_records = n;
    return new_buffer;
}

