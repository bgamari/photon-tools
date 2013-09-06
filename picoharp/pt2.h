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

#pragma once
#include <stdint.h>
#include <stdexcept>
#include <istream>
#include <vector>
#include "picoharp.h"

#define PT2_TIME_UNIT 4E-12 // 4 ps
#define PT2_WRAPAROUND_TIME 0x0c8f0000 

/* T2 mode event record */
struct pt2_record { 
	uint8_t channel;
	uint64_t time;
	bool special;
	uint8_t markers;
};

class pt2_file : public picoharp_file {
private:
	unsigned int nsync;
public:
	pt2_file(picoharp_file& pf) : picoharp_file(pf), nsync(0) {
		if (binary_hdr.meas_mode != PT2_MEASMODE_T2)
			throw std::runtime_error("Unsupported measurement mode");
	}
	pt2_file(std::istream& is) : pt2_file(is) {}
	pt2_record read_record();
	std::vector<pt2_record> read_all_records();
};

// A convenient helper
uint64_t *get_pt2_timestamps(std::istream& is,
                             unsigned int channel,
                             unsigned int *n_records);

