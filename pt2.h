#define PT2_DISPCURVES 8
#define PT2_TIME_UNIT 4E-12 // 4 ps
#define PT2_WRAPAROUND_TIME 210698240 
#define PT2_MEASMODE_T2 2 


class pt2_file {
private:
	std::istream& is;
	uint32_t overflow_time;
	uint32_t records_read;

public:
	pt2_text_hdr text_hdr;
	pt2_binary_hdr binary_hdr;
	pt2_board_hdr board_hdr;
	pt2_tttr_hdr tttr_hdr;
	std::vector<pt2_record>* records;


	pt2_file(std::istream& is) : is(is), overflow_time(0), records_read(0), records(NULL) {
		read_headers();
	}

	pt2_record read_record() {
		uint32_t rec;
		pt2_record r;

		if (records_read >= tttr_hdr.n_records)
			throw new ArgumentException("Read too many records");
		records_read++;

		is.read(&rec, 4);
		r.time = 0x0fffffff & rec;
		r.channel = 0xf0000000 & rec;
		
		if (r.channel == 0xf) {
			r.special = true;
			r.markers = r.time & 0xf;
			if (markers == 0) // overflow record
				overflow_time += PT2_WRAPAROUND_TIME;
		}

		r.time += overflow_time;

		return r;
	}

	void read_all_records(istream& is) {
		int n = tttr_hdr.n_records;
		records = new std::vector<pt2_record>(n);
		for (int i=0; i < n; i++)
			records.push_back(read_record());
	}

private:
	void read_headers() {
		is.read(&text_hdr, sizeof(pt2_text_hdr));
		check_text_header();
		is.read(&binary_hdr, sizeof(pt2_binary_hdr));
		check_binary_header();
		is.read(&board_hdr, sizeof(pt2_board_hdr));
		is.read(&tttr_hdr, sizeof(pt2_tttr_hdr));
		is.ignore(4*tttr_hdr.imaging_hdr_sz);
		read_all_records(is);
	}


	void check_text_header() {
		if (strcmp(text_hdr.ident, "PicoHarp 300"))
			throw new ArgumentException("File identifier not found");

		if (strncmp(text_hdr.format_version, "2.0", 3))
			throw new ArgumentException("Unsupported file format version");
	}

	void check_binary_header() {
		if (binary_hdr.meas_mode != PT2_MEASMODE_T2)
			throw new ArgumentException("Unsupported measurement mode");
	}
};


#pragma pack(4)

/* ASCII portion of PT2 header */
struct pt2_text_hdr {
	char ident[16];			//"PicoHarp 300"
	char format_version[6];		//file format version
	char creator_name[18];		//name of creating software
	char creator_version[12];	//version of creating software
	char file_time[18];
	char crlf[2];
	char comment[256];
};

/* Binary portion of PT2 header */
struct pt2_binary_hdr {
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
	} disp_curves[PT2_DISPCURVES];	

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
struct pt2_board_hdr {		
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
struct pt2_tttr_hdr {
	uint32_t ext_devices;
	uint32_t reserved1;
	uint32_t reserved2;			
	uint32_t counter_rate[2];
	uint32_t stop_after;
	uint32_t stop_reason;
	uint32_t n_records;
	uint32_t imaging_header_sz;
};


/* T2 mode event record */
struct pt2_record { 
	uint8_t channel;
	uint64_t time;
	bool special;
	uint8_t markers;
};

