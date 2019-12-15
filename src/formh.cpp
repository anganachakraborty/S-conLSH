
#include "formh.h"
#include "desc.h"


#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>


using namespace std;

char *const short_options = "h:K:L:z:l:w:m:x:y:q:r:a:n";
struct option long_options[] = {
    
    { "help",     0,   NULL,    'h'   },  
    {"K",1,NULL,'K'},
    {"L",1,NULL,'L'},
    {"zero",1,NULL,'z'},
    {"lambda",1,NULL,'l'},
    { "window-hits",     1,   NULL,    'w'   },
    { "candidates",     1,   NULL,    'm' },
    
    { "match",     1,   NULL,    'x'   },
    { "mismatch",		1,NULL,	'y'},
    
    {"gap-open", 1,  NULL, 'q'},
    {"gap-extension", 1,  NULL,'r'},
    {"align", 1,  NULL,'a'},
   // {"threads", 1, NULL, 't'},
    //{"local-kmer", 1, NULL, 'l'},
    //{"help",    0,  NULL,'h'},
    { 0,     0,   0,    0   }
};

Form::Form(opts *opt)
{
 // no of 1's max 16
    opt->len_sed = 32;
    
    opt->K=2;
    opt->L=1;
    opt->zero=5;
    opt->lambda=7; // context size 2*lambda+1
    //opt->len_sed = 32;
    opt->canN = 800;
    opt->hit_limit = 1000;
    opt->thread = 1;
    opt->len_limit = 100000;
    opt->gapopen = 2;
    opt->gapextend = 1;
    opt->mismatch = 5;
    opt->match = 2;
    opt->algn=0;
   
}

int Form::usage()
{

    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   conLSH-indexer\n"); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n\n", CONTACT); 
    fprintf(stderr, "Usage:     conLSH-indexer [Options] <HashIndexDir> <Reference> <ReadFile>\n\n"); 

    fprintf(stderr, "<HashIndexDir>         The directory storing conLSH index\n");
    fprintf(stderr, "<Reference>            Sequence of reference genome, in FASTA format\n\n");
    fprintf(stderr, "<ReadFile>             Reads file, in FASTQ/FASTA format\n");
    fprintf(stderr, "           -K, --K                <int>           Concatenation factor of locality sensitive hashing\n"); 
    fprintf(stderr, "           -L, --L                <int>           No. of hash tables in conLSH framework\n");  
    fprintf(stderr, "           --lambda               <int>           The context factor\n"); 
    fprintf(stderr, "           --zero                 <int>           The number of don't cares or zeros in the S-conLSH pattern\n"); 
    fprintf(stderr, "           -w, --window-hits      <int>           the max allowed number of windows hitting by a k-mer [1000]\n"); 
    fprintf(stderr, "           -m, --candidates       <int>           the number of candidates for extension [400]\n"); 
   	fprintf(stderr, "           -x, --match            <int>           score of match for the alignments in extension phase [2]\n");
   	fprintf(stderr, "           -y, --mismatch         <int>           mismatch penalty for the alignments in extension phase [5]\n");
   	fprintf(stderr, "           -q, --gap-open         <int>           gap open penalty for the alignments in extension phase [2]\n");
   	fprintf(stderr, "           -r, --gap-extension    <int>           gap extension penalty for the alignments in extension phase [1]\n");
   	fprintf(stderr, "           -a, --align		  <int>           Value=1, outputs alignment in SAM format[Default=0, Alignment Free PAF format output]\n");
   	//fprintf(stderr, "           -l, --local-kmer       <int>           the minimum length of the local matches used for SDP [11]\n");
    
    fprintf(stderr, "           -h, --help                             help\n");
     
    fprintf(stderr, "\n"); 
    return 0;
}

int Form::opt_parse(int argc, char *argv[], opts* opt)
{
	int c; 
    int option_index=0;
    if (argc == 1) return usage();
    while((c = getopt_long(argc, argv, short_options, long_options, &option_index))>=0){
        switch(c){
            case 'h':
                return usage();
                break;
            case 'K':
                opt->K = atoi(optarg);
                break;
            case 'L':
                opt->L = atoi(optarg);
                break;
            case 'z':
                opt->zero = atoi(optarg);
                break;
            case 'l':
                opt->lambda = 2*atoi(optarg)+1;
                break;
            case 'm':
                opt->canN = atoi(optarg);
                break;
            case 'w':
                opt->hit_limit = atoi(optarg);
                break;
            
            case 'q':
                opt->gapopen = atoi(optarg);
                break;
            case 'r':
                opt->gapextend = atoi(optarg);
                break;
            case 'x':
                opt->match = atoi(optarg);
                break;
            case 'y':
                opt->mismatch = atoi(optarg);
                break;
            case 'a':
                opt->algn = atoi(optarg);
                break;
            default:
                fprintf(stderr,"not proper parameters\n");
                return usage();
              
        }
    
    }
    if(optind + 3 != argc){
         fprintf(stderr, "[opt_parse]: index directory, read file and reference file can't be omited!\n"); 
        return 0; 
    }
  
 	strncpy(opt->hashdir, argv[optind++],sizeof(opt->hashdir));

 	strncpy(opt->refpath,argv[optind++],sizeof(opt->refpath));
    strncpy(opt->readpath,argv[optind++],sizeof(opt->readpath));
    return 1;  

}
