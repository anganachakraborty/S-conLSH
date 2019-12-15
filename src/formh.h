#include "desc.h"
#include <stdint.h>

#define PATH_LEN 1024
#define LEN_BASES 1024

#define LEN 100

typedef struct options {
	uint32_t 	len_sed;
	char 		refpath[PATH_LEN];
	char 		hashdir[PATH_LEN];
	char 		readpath[PATH_LEN];
	//uint32_t		w;
	uint32_t 		K;
	uint32_t 		L;
	uint32_t 		lambda; // context factor
	uint32_t		zero;
	uint32_t 	canN;
	uint16_t 	hit_limit;
	uint32_t 	len_limit;
        int 		thread;
	int 		argc;
	char 		**argv;
	int 		gapopen;
	int		gapextend;
	int		match;
	int 		mismatch;
	int		algn;
	

}opts;


class Form {
	
	opts *opt;
public:
	Form(opts *opt);
	int usage();
	int opt_parse(int argc, char *argv[], opts* opt);
};
