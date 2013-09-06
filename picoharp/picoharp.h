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
#include <istream>

#define PT_DISPCURVES 8

enum measurement_mode {
        PT2_MEASMODE_T2 = 2,
        PT2_MEASMODE_T3 = 3,
};

/* ASCII portion of PT2/3 header */
struct pt_text_hdr {
	char ident[16];			//"PicoHarp 300"
	char format_version[6];		//file format version
	char creator_name[18];		//name of creating software
	char creator_version[12];	//version of creating software
	char file_time[18];
	char crlf[2];
	char comment[256];
};

/* Binary portion of PT2/3 header */
struct pt_binary_hdr {
	uint32_t n_curves;
	uint32_t bits_per_record;
	uint32_t routing_channels;
	uint32_t n_boards;
	uint32_t active_curve;
	uint32_t meas_mode;
	uint32_t sub_mode;
	uint32_t range;
	uint32_t offset;
	uint32_t acq_time;			// in ms
	uint32_t stop_at;
	uint32_t stop_on_overflow;
	uint32_t restart;
	uint32_t disp_lin_log;
	uint32_t disp_time_from;		// 1ns steps
	uint32_t disp_time_to;
	uint32_t disp_counts_from;
	uint32_t disp_counts_to;

	struct {
		uint32_t map_to;
		uint32_t show;
	} disp_curves[PT_DISPCURVES];	

	struct {
		float start;
		float step;
		float end;
	} params[3];

	uint32_t repeat_mode;
	uint32_t repeats_per_curve;
	uint32_t repeat_time;
	uint32_t repeat_wait_time;
	char script_name[20];
};

/* Board header */
struct pt_board_hdr {		
	char hardware_ident[16]; 
	char hardware_version[8]; 
	uint32_t hardware_serial; 
	uint32_t sync_divider;
	struct {
		uint32_t zero_cross;
		uint32_t level;
	} cfds[2];
	float resolution;

	// below is new in format version 2.0
	uint32_t router_model_code;
	uint32_t router_enabled;
	struct {
		uint32_t input_type;
		uint32_t input_level;
		uint32_t input_edge;
		uint32_t cfd_present;
		uint32_t cfd_level;
		uint32_t cfd_zero_cross;
	} router_channels[4];
};

/* TTTR mode header */
struct pt_tttr_hdr {
	uint32_t ext_devices;
	uint32_t reserved1;
	uint32_t reserved2;			
	uint32_t counter_rate[2];
	uint32_t stop_after;
	uint32_t stop_reason;
	uint32_t n_records;
	uint32_t imaging_hdr_sz;
};

class picoharp_file {
protected:
	std::istream& is;
	uint64_t overflow_time;
	uint32_t records_read;

public:
	pt_text_hdr text_hdr;
	pt_binary_hdr binary_hdr;
	pt_board_hdr board_hdr;
	pt_tttr_hdr tttr_hdr;

protected:
	picoharp_file(std::istream& is) : is(is), overflow_time(0), records_read(0) {
		read_headers();
	}

private:
	void read_headers();

public: 
        static uint64_t* get_timestamps(std::istream& is,
                                        unsigned int channel,
                                        unsigned int* n_records);
};
