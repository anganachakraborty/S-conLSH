#include <stdio.h>
#include "formh.h"
#include "whash.h"
#include "readfl.h"
#include <iostream>
#include <time.h>
#include <cstdlib>
#include <sys/resource.h>
#include "aligner.h"

using namespace std;


const std::string getCurrentDateTime(time_t *now) {
    *now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(now);
    strftime(buf, sizeof(buf), "[INFO] %Y-%m-%dT%X", &tstruct);

    return buf;
}

int main(int argc, char *argv[])
{
	
	time_t time1, time2, time3, time4;
	
	
	opts *opt = new opts;
	Form fm(opt);
	if (fm.opt_parse(argc,argv,opt)!=1)
		exit(1);

	
	time(&time1);
	//fprintf(stderr,"%s conLSH-indexer started\n",getCurrentDateTime(&time1).c_str());
	
	uint32_t len_genome = 0;
	int 		countChr;
	char **ChrName;
	uint32_t 	Start_pos[100];
	read_file rd;
	//char 	*genome = rdfl.read_ref(opt->refpath,&len_genome);
	ChrName = new char *[LEN];
	for (int i=0;i<100;++i) { ChrName[i] = new char[LEN];}
	//int 	 countChr;
	fprintf(stderr,"Reference genome indexing.....\n");
	char *genome = rd.read_ref(opt->refpath,&len_genome,Start_pos,ChrName,&countChr);
	if (NULL == genome) {

		fprintf(stderr,"Fail to load reference genome, now exit");
		exit(1);
	}
	else
	  fprintf(stderr,"Genome Length = %d\nNumber of chromosomes = %d\n\n",len_genome,countChr);
	char *genome_e = genome + len_genome - 1;

	if (NULL == genome) {

		fprintf(stderr,"Fail to load Genome reference, now exit");
		exit(1);
	}

	Hash hashh;
	
	hashh.write_hashfile(opt->hashdir,genome,len_genome,opt->len_sed,opt->K,opt->L,opt->lambda, opt->zero);
	 time(&time2);
   struct rusage r_usage;
      
	//fprintf(stderr,"%d\t%d\t%d\t%d\t%.2f",opt->K,opt->lambda,opt->L,opt->zero,difftime(time2, time1));
	fprintf(stderr,"Parameters used...\n K=%d\n Context size (2*lambda+1) = %d\n L=%d\n zero=%d\n\nIndexing time %.0f sec\n\n\n",opt->K,opt->lambda,opt->L,opt->zero,difftime(time2, time1));
    
	time(&time3);
	fprintf(stderr,"Aligning...\n ");
		
	Aligner alig(opt->readpath,opt->thread,argc,argv,opt->gapopen,opt->gapextend,opt->match,opt->mismatch,opt->canN,len_genome,countChr,ChrName,Start_pos,genome,genome_e, hashh.pat, opt->L, hashh.h_index, hashh.hv, opt->algn);
	alig.Runtask();	

	
	time(&time4);
	
	 
    getrusage(RUSAGE_SELF,&r_usage);
  
	fprintf(stderr,"Parameters used...\n Match Score=%d\n Mismatch penalty=%d\n Gap_open penalty=%d\n Gap_extend penalty=%d\n\nAlignment time %.0f sec\nPeak memory foorprint %ld MB\n",opt->match,opt->mismatch,opt->gapopen,opt->gapextend,difftime(time4, time3),r_usage.ru_maxrss/1000);
	
	if ( NULL != opt ) delete opt;
	
	if (NULL != ChrName) {
		for (int i=0;i<100;++i) { delete ChrName[i]; }
		delete ChrName;
	}
	return 0;

}
