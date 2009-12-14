#include "pt2.h"
#include <stdexcept>
#include <cstring>

pt2_record pt2_file::read_record() {
	uint32_t rec;
	pt2_record r;

	records_read++;
	if (records_read > tttr_hdr.n_records)
		throw new std::runtime_error("Read too many records");

	is.read((char*) &rec, 4);
	r.time = 0x0fffffff & rec;
	r.channel = 0xf0000000 & rec;
	
	if (r.channel == 0xf) {
		r.special = true;
		r.markers = r.time & 0xf;
		if (r.markers == 0) // overflow record
			overflow_time += PT2_WRAPAROUND_TIME;
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
		throw new std::runtime_error("File identifier not found");
	if (strncmp(text_hdr.format_version, "2.0", 3))
		throw new std::runtime_error("Unsupported file format version");

	is.read((char*) &binary_hdr, sizeof(pt2_binary_hdr));
	if (binary_hdr.meas_mode != PT2_MEASMODE_T2)
		throw new std::runtime_error("Unsupported measurement mode");

	is.read((char*) &board_hdr, sizeof(pt2_board_hdr));
	is.read((char*) &tttr_hdr, sizeof(pt2_tttr_hdr));
	is.ignore(4*tttr_hdr.imaging_hdr_sz);
}

