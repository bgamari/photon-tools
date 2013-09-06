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

#include "picoharp.h"
#include "pt2.h"
#include "pt3.h"
#include <cstring>

void picoharp_file::read_headers()
{
	is.read((char*) &text_hdr, sizeof(pt_text_hdr));
	if (strcmp(text_hdr.ident, "PicoHarp 300"))
		throw std::runtime_error("File identifier not found");
	if (strncmp(text_hdr.format_version, "2.0", 3))
		throw std::runtime_error("Unsupported file format version");

	is.read((char*) &binary_hdr, sizeof(pt_binary_hdr));
	is.read((char*) &board_hdr, sizeof(pt_board_hdr));
	is.read((char*) &tttr_hdr, sizeof(pt_tttr_hdr));
	is.ignore(4*tttr_hdr.imaging_hdr_sz);
}

uint64_t* picoharp_file::get_timestamps(unsigned int channel,
                                        unsigned int* n_records)
{
        switch (binary_hdr.meas_mode) {
        case PT2_MEASMODE_T2:
                return pt2_file(*this).get_timestamps(channel, n_records);
        case PT2_MEASMODE_T3:
                return pt3_file(*this).get_timestamps(channel, n_records);
        default:
                return NULL;
        }
}
