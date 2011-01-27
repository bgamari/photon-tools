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

#include <boost/format.hpp>
#include <iostream>
#include <fstream>
#include "pt2.h"

// channel=0xf dumps all channels
void dump(std::istream& is, std::ostream& os, uint8_t channel=0xf) {
        pt2_file pt2(is);
	unsigned int n_rec = pt2.tttr_hdr.n_records;

	for (int i=0; i < n_rec; i++) {
		pt2_record rec = pt2.read_record();
		if (!rec.special && (channel == 0xf || channel == rec.channel)) {
			uint64_t time = rec.time;
			int count = 1;
			os.write((char*) &time, sizeof(uint64_t));
		}
	}
}

int main(int argc, char** argv) {
        std::string out_name;
        if (argc > 1) {
                std::string name = argv[1];
                std::ifstream is(name);
                int chan;

                if (argc > 2) {
                        chan = atoi(argv[2]);
                        out_name = (boost::format("%s.ch%d.times") % name % chan).str();
                } else {
                        chan = 0xf;
                        out_name = (boost::format("%s.times") % name).str();
                }

                std::ofstream os(out_name);
                dump(is, os, chan);

                std::ofstream pos(out_name + ".meta");
                pos << (boost::format("{\n"
                                        "\"clockrate\": %f\n"
                                        "\"date\": \"%s\"\n"
                                      "}") % (1.0/PT2_TIME_UNIT));
        } else
                dump(std::cin, std::cout);

	return 0;
}

