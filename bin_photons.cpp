/* timetag-tools - Tools for UMass FPGA timetagger
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


#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <boost/program_options.hpp>

/*
 *
 * Temporally bins a photon stream
 *
 * Usage:
 *   bin_photons [BIN_LENGTH]
 *
 * Where BIN_LENGTH is the length of each bin in counter units.
 *
 * Input:
 *   A binary photon stream
 *
 * Output:
 *   A binary stream of bin_records
 *
 * Notes:
 *   We handle wrap-around here by simply keeping all times as 64-bit and
 *   assume that we wrap-around at most once.  With 1 nanosecond clock units,
 *   this gives us 500 years of acquisition time.
 *
 */

void write_bin(uint64_t time, uint16_t counts) {
        using std::cout;
        cout.write((char*) &time, 8);
        cout.write((char*) &counts, 2);
}

int main(int argc, char** argv) {
        namespace po = boost::program_options;
        using std::cin;
        using std::cout;

        po::options_description opts;
        opts.add_options()
                ("help", "produce help message")
                ("no-buffer", "disable buffering of output (for realtime use)")
                ("no-zeros", "disable output of bins with no photons")
                ("bin-length", po::value<uint64_t>(), "length of bins");
        po::positional_options_description popts;
        popts.add("bin-length", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv)
                        .options(opts).positional(popts).run(), vm);
        po::notify(vm);

        if (vm.count("help") || !vm.count("bin-length")) {
                cout << opts << "\n";
                return 1;
        }

        bool no_zeros = vm.count("no-zeros");
        uint64_t bin_length = vm["bin-length"].as<uint64_t>();
        if (vm.count("no-buffer"))
                setvbuf(stdout, NULL, _IONBF, 0);

        // Read first photon for start time
        uint16_t bin_count = 1;
        uint64_t bin_start;
        cin.read((char*) &bin_start, 8);

        while (cin.good()) {
                uint64_t time;
                cin.read((char*) &time, 8);
                if (time >= bin_start+bin_length) {
                        uint64_t new_bin_start = (time / bin_length) * bin_length;

                        // First print photons in last bin
                        write_bin(bin_start, bin_count);

                        // Then print bins with no photons since last bin
                        if (!no_zeros) {
                                for (uint64_t t=bin_start+bin_length; t < new_bin_start; t += bin_length)
                                        write_bin(t, 0);
                        }
                        
                        // Then start new bin
                        bin_count = 0;
                        bin_start = new_bin_start;
                }
                bin_count++;
        }

        if (!cin.eof()) return 2;
        return 0;
}

