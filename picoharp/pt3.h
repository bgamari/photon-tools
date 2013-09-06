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

#pragma once
#include <stdint.h>
#include <stdexcept>
#include <istream>
#include <vector>
#include "picoharp.h"

#define PT3_WRAPAROUND_TIME 0x00010000 
:
/* T3 mode event record */
struct pt3_record { 
	bool is_special;
	uint16_t numsync;
	uint8_t channel;
	union {
		struct {
			uint64_t time;
		} normal;
		struct {
			uint8_t markers;
		} special;
	};
};

class pt3_file : public picoharp_file {
	unsigned int nsync;
	unsigned int sync_period; // in cycles
public:
	pt3_file(std::istream& is) : pt3_file(is) { }
	pt3_file(picoharp_file& pf) :
		picoharp_file(pf),
		nsync(0),
		sync_period(1e-9 * tttr_hdr.counter_rate[0] / board_hdr.resolution)
	{
		if (binary_hdr.meas_mode != PT2_MEASMODE_T3)
			throw std::runtime_error("Unsupported measurement mode");
	}
	pt3_record read_record();
	std::vector<pt3_record> read_all_records();
};

// A convenient helper
uint64_t *get_pt3_timestamps(const char *filename,
                             unsigned int channel,
                             unsigned int *n_records);

