
#ifndef ALIGNER_H_
#define ALIGNER_H_
#define LEN_BASES 1024

#define LEN 100

#define MAX_SN 6

//#include "graph.h"
#include "whash.h"
#include "kseq.h"
#include "ksw.h"
#include "desc.h"
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)

#include <iostream>
#include <cstdio>
#include <cstring>
#include <queue>
using namespace std;

const uint8_t seq_nt4_tablet[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};




typedef struct sam_rec{
	uint8_t 	chrIndex;
	uint32_t 	pos;
	//char 		_cigar[LEN_LIMIT<<1];
	//string		headCigar;
	//string 		bodyCigar;
	//string 		tailCigar;
	string		cigar;
	int 		score;

	//int 		headScore;
	//int 		bodyScore;
	//int 		tailScore;
	//int 		score;
	uint16_t 	flag;
	int16_t		MAQ;
	uint32_t        chrLen;
	struct sam_rec &	operator=(const struct sam_rec &r)  
	{ 
		chrIndex = r.chrIndex; 
		pos = r.pos;
		cigar = r.cigar;
		score = r.score;
		flag = r.flag;
		MAQ = r.MAQ;
		return *this;
	}
}Sam_Rec;


typedef struct svsam_rec {
	uint8_t 	chrIndex;
	uint32_t 	pos;
	//char 		_cigar[LEN_LIMIT<<1];
	string		cigar;
	int 		score;
	uint16_t 	flag;
	
	uint32_t 	ref_end;
	uint16_t 	read_end;
	uint16_t	read_start;
	uint32_t 	ref_start;
	uint16_t	lclip;
	uint16_t 	rclip;
	int16_t 	MAQ;
	struct svsam_rec &	operator=(const struct svsam_rec &r)  
	{ 
		chrIndex = r.chrIndex; 
		pos = r.pos;
		// headCigar = r.headCigar;
		// bodyCigar = r.bodyCigar;
		// tailCigar = r.tailCigar;
		cigar = r.cigar;
		
		score = r.score;
		flag = r.flag;
		ref_start = r.ref_start;
		rclip = r.rclip;
		lclip = r.lclip;
		MAQ = r.MAQ;
		return *this;
	}
}SvSam_Rec;


typedef struct bucket2
{
	uint16_t	*hit_times;
	uint32_t 	*seq_num;
	bool		*isrc;
	
}bkt2;

typedef struct {
	char	*qual;
	int 	qual_len;

	char   	*seq;
	char 	*rseq;
	int 	seq_len;

	char 	*name;
	int 	name_len;
}rhat_seq;

typedef struct {
	uint32_t		*sed_rec;
	uint16_t 		*sed_hit_times;
	uint16_t		*unused_bkt;

	//RHashtable 		*rhashtab;
	//RHashtable 		*rrhashtab;

}aux_var;

class Aligner {
	bkt2 		preserved;// store the target sites
	//uint32_t 	pos[20];//write	
	//uint8_t 	chrIndex[20];//write
	int 		countChr;//readable  
	char 		*genome;//readable
	char 		*genome_e;
	uint32_t 	len_genome;//readable
	//Hashtab 	*hashtab;//readable
	char 		**ChrName;//readable
	uint32_t 	*Start_pos;//readable
	//opts 		*opt;//readable
	int8_t 		mat[25];
         int 		gapopen;
	int		gapextend;
	int		match;
	int 		mismatch;
	uint32_t 	canN;
        char 		*readpath;
        int 		thread;
	int 		argc;
	char 		**argv;
	uint32_t 	*pat;
	uint32_t 	L;
	 H_index * h_index;
	uint32_t *hv;
	int	algn;
	//uint32_t 	*hhv;
	//uint32_t 	*h_index;
public:
	Aligner(char *rpath, int thrd, int ac, char **av, int gopen,int gextend,int mtch,int mismtch, uint32_t cN, uint32_t l_genome, int cChr, char **CName, uint32_t *s_pos, char *g, char *g_e, uint32_t *p, uint32_t LL, H_index * h_indx,uint32_t *hhv, int aa);
	~Aligner();
	void Runtask();
	
	int 			applyNonSV(kseq_t *trunk, Sam_Rec *_sams, uint32_t L, uint32_t *pat, uint32_t len_genome, int align, int tr);

	int 			applySV(kseq_t *trunk, Sam_Rec *_svsams, uint32_t L, uint32_t *pat, uint32_t len_genome, int align, int tr);
private:
	uint32_t 		*sed_rec;//len_bases -14] ;//use malloc();
	uint16_t 		*sed_hit_times;//use malloc
	uint16_t 		*unused_bkt;
	
	//uint32_t 		sv_interval;
	

	//int 			conductAlign(kseq_t *trunk, char *read, char *rcRead, int lenRead, RHashtable *_rhashtab, RHashtable *_rrhashtab
	//,int c, SvSam_Rec *_svsamsp); 
	int 			align_NonSV(kseq_t *trunk, char *useread, Sam_Rec *_sams);
	int 			align_SV(kseq_t *trunk, int offset, char *useread,Sam_Rec *_sams);

	int 			alignFree_NonSV(kseq_t *trunk, char *useread ,Sam_Rec *_sams);
	int 			alignFree_SV(kseq_t *trunk, int offset, char *useread, Sam_Rec *_sams);
	//int 			conductAlign(kseq_t *trunk,int c, Sam_Rec *_sams); 
	
	//void 			proCans(char *read, uint32_t lenRead, bool isRC,  std::priority_queue<bkt2> &cansHeap, 
	//uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt) ;

	//uint32_t 		gen_sed(char *bases, uint32_t len_bases, uint32_t *_sed_rec);
	
	//void 			sort_sed(uint32_t *_sed_rec,uint32_t usedseed);

	//void  			stat_sed(uint32_t *_sed_rec, uint16_t *_sed_hit_times, uint16_t *_unused_bkt, uint32_t usedseed);
	
	int 			revComRead(char *read, char *rcRead, int len_read);
	
	bool 			isQualifiedRead(char *read,int len);

	
	//sv
	//int 			rhat_seq_read(kstream_t *_fp, kseq_t *_seqs, int n_needed);
	//int 			qualified(SvSam_Rec *key,SvSam_Rec *set, int start, int end);
	//int 			connect(SvSam_Rec *rec1, SvSam_Rec *rec2, kseq_t *trunk);
	//int 			produceSAM(SvSam_Rec *_svsams , int countbulks,int *sam4bulk, kseq_t *trunk, uint32_t *len);
	int 			mergeBoundaryCigar(string & cigar1, string & cigar2, uint32_t *cigar, int n_cigar, const char *correspondTable);
	//output
	//int 			outputSam(kseq_t *_seqs, Sam_Rec *_sams, SvSam_Rec **_svsams, uint16_t *_sam_details, int _n_seqs);
};
typedef struct aux
{
	int 		tid;
	Aligner 	*aln;
	int 		n_seqs;
	aux_var 	com_var;
	kseq_t 		*seqs;
	Sam_Rec		*sams;
	SvSam_Rec	**svsams;
	//opts 		*opt;
	uint16_t 	*sam_details;// H:6 samamount L:2 isSAM, bits set by limit of read length might be changed in future
}thread_aux;
//thread_aux *thread_initiate(int n_thread, RHashtable **rhashtab, RHashtable **rrhashtab, uint32_t *sed_rec, uint16_t *sed_hit_times, 
	//uint16_t *unused_bkt, opts *_opt, Aligner *aln) ;
//static void 	*thread_worker(void *data);
extern 	uint8_t rev[];
int 	transIntoDec(uint8_t *transstr,char *str,int length);
extern const uint8_t seq_nt4_tablet[];
#endif
