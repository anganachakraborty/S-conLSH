
#include "aligner.h" 

#include "conLSH.h"
#include "whash.h"
#include <pthread.h> 
//#include "btree.h"


#define SELECT_NUM 10

#define N_NEEDED 5000

#define FORELEN  5




uint8_t rev[128]={
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
	4,84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,
	4,4,4,4,65,4,4,4,4,4,4,4,4,4,4,4,4,
	84,4,71,4,4,4,67,4,4,4,4,4,4,67,4,4,
	4,4,4,65,4,4,4,4,4,4,4,4,4,4,4
};




int read_seq;
pthread_rwlock_t rwlock;
//BTree t(3); 

void revstr(uint8_t *revstring, char *string, int len)
{
	for (int i=0;i<len;++i) {
		revstring[i] = seq_nt4_tablet[string[len - 1 - i]];
	}

}
int compare_bkt2(const void *p,const void *q)
{
	const bkt2 *f = (bkt2 *)p;
	const bkt2 *t = (bkt2 *)q;
	return (int)(f->seq_num - t->seq_num);
}

int compare_sam(const void *p, const void *q, void *arg)
{
	int f = *(int *)p;
	int h = *(int *)q;
	Sam_Rec *s = (Sam_Rec *)arg;
 	return s[h].score - s[f].score;
}
//void copy_array_sorted(uint32_t *a, int n1, uint32_t *b, int n2)
//{
	////INPUT: a and b are two sorted arrays
	////OUTPUT: a sorted array of (n1+n2) elements
	
	
	
	//int p1, p2, k;
	//p1=n1-1;
	//p2=n2-1;
	//k=(n1+n2-1);
	
	//while(p1>=0 && p2>=0)
	//{
		//if(a[p1]>b[p2])
		//{
			//a[k]=a[p1];
			//p1--;
			//k--;
		//}
		//else 
		//{
			//a[k]=b[p2];
			//p2--;
			//k--;
		//}
		
	//}
	//while(p1>=0)
	//{
		//a[k]=a[p1];
			//p1--;
			//k--;
	//}
	//while(p2>=0)
	//{
		//a[k]=b[p2];
			//p2--;
			//k--;
	//}
	
	
	
//}

Aligner::Aligner(char *rpath, int thrd, int ac, char **av, int gopen,int gextend,int mtch,int mismtch, uint32_t cN,uint32_t l_genome,int cChr,char **CName,uint32_t *s_pos, char *g, char *g_e, uint32_t *p, uint32_t LL, H_index * h_indx,uint32_t *hhv, int aa)
{
	//opt = _opt;
	len_genome=l_genome;
	countChr=cChr;
	ChrName=CName;
	Start_pos=s_pos;
	genome=g;
	genome_e=g_e;
	readpath=rpath;
	thread=thrd;
	argc=ac;
	argv=av;
	gapopen=gopen;
	gapextend=gextend;
	match=mtch;
	mismatch=mismtch;
	canN=cN;
	pat=p;
	L=LL;
	h_index=h_indx;
	hv=hhv;
	algn=aa;
	//cout<<"in aligner"<<hv[80515860]<<endl;
	//sed_rec = new uint32_t[opt->thread * LEN_BASES];
	//if(NULL == sed_rec) {
		//fprintf(stderr, "Failed when applying for new space! now exit");
		//exit(1);
	//}

	//sed_hit_times = new uint16_t[opt->thread * LEN_BASES];
	//if (NULL == sed_hit_times) {
		//fprintf(stderr, "Failed when applying for new space! now exit");
		//exit(1);
	//}

	//unused_bkt = new uint16_t[opt->thread * LEN_BASES];
	//if (NULL == unused_bkt) {
		//fprintf(stderr, "Failed when applying for new space! now exit");
		//exit(1);
	//}

	// RCRead = new char[opt->thread * opt->len_limit];
	// if (NULL == RCRead) {
	// 	fprintf(stderr, "Failed when applying for new space! now exit");
	// 	exit(1);
	// }

	//initiate ksw parameters
	int k;
	for (int i = k = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
			mat[k++] = i == j? match : -mismatch;
		mat[k++] = 0; // ambiguous base
	}
	for (int j = 0; j < 5; ++j) mat[k++] = 0;

	
}

Aligner::~Aligner()
{
	//if (NULL != sed_rec)
		//delete[] sed_rec;
	//if (NULL != sed_hit_times)
		//delete[] sed_hit_times;
	//if (NULL != unused_bkt)
		//delete[] unused_bkt;
	//// if (NULL != RCRead)
	// 	delete[] RCRead;
}

//status: checked
//problem: none
int Aligner::revComRead(char *read, char *rcRead, int len_read)
{
	int i;
	//uint32_t offset = groupNum * opt->len_limit;
	for(i=0;i<len_read;i++)
		rcRead[i] = rev[read[len_read - 1 - i]];
	rcRead[i] = '\0';
	return 0;
}


//status: checked
//problem: none
int compar(const void *p,const void *q)
{
	const uint32_t *t = (uint32_t *)p;
	const uint32_t *f = (uint32_t *)q;
	if(*t>*f)
		return 1;
	else
		if (*t<*f)
			return -1;
		else
			return 0;
}





uint32_t fig_pos(uint32_t pos,uint32_t offset)
{
	uint32_t bkt_num = pos;
	if (bkt_num <= offset) {
		return 0;
	}

	uint32_t left_start =  bkt_num - offset ;
	return left_start;
}

int 	transIntoDec(uint8_t *transtr,char *str,int length)
{
	for (int i=0;i<length;++i) {
		transtr[i] = seq_nt4_tablet[str[i]];
	}
	return 0;
}
int Aligner::mergeBoundaryCigar(string & cigar1, string & cigar2, uint32_t *cigar, int n_cigar, const char *correspondTable) 
{
	uint32_t len_cigar1 = cigar1.size();
	uint32_t len_cigar2 = cigar2.size();
	
	uint32_t digit_end_cigar2,digit_start_cigar1;

	for (digit_end_cigar2 = 0; digit_end_cigar2 < len_cigar2; ++digit_end_cigar2) {
		if (isalpha(cigar2[digit_end_cigar2]))
			break;
	}

	for (digit_start_cigar1 = len_cigar1 -2; digit_start_cigar1 >0;--digit_start_cigar1) {
		if (isalpha(cigar1[digit_start_cigar1])) 
			break;
	}

	uint32_t len_digit = len_cigar1 - digit_start_cigar1 - 1; 
	uint32_t countAlpha_cigar1 = strtoul(cigar1.substr(digit_start_cigar1 + 1, len_digit).c_str(),0,10);
	uint32_t countAlpha_cigar2 = strtoul(cigar2.substr(0, digit_end_cigar2).c_str(),0,10);

	
	if (correspondTable[cigar[0]&0xf] ==cigar1[len_cigar1-1]) {
		cigar[0] += (countAlpha_cigar1 << 4);
		cigar1.erase(digit_start_cigar1+1,len_digit + 1);
	}

	if ( correspondTable[cigar[n_cigar-1]&0xf] == cigar2[digit_end_cigar2]) {
		cigar[n_cigar-1] += (countAlpha_cigar2 << 4) ;

		cigar2.erase( 0, digit_end_cigar2 + 1);
	}

	return 0;
}










/*thread_aux *thread_initiate(int n_thread, RHashtable **rhashtab, RHashtable **rrhashtab, uint32_t *sed_rec, uint16_t *sed_hit_times,
	uint16_t *unused_bkt, opts *_opt, Aligner *aln)
{
	thread_aux *aux = new thread_aux[n_thread];

	for (int i= 0; i<n_thread; ++i) {
		aux[i].aln = aln;
		aux[i].opt = _opt;
		aux[i].com_var.sed_hit_times = sed_hit_times + i * LEN_BASES;
		aux[i].com_var.unused_bkt = unused_bkt + i * LEN_BASES;
		aux[i].com_var.sed_rec = sed_rec + i * LEN_BASES;
		aux[i].com_var.rhashtab = rhashtab[i];
		aux[i].com_var.rrhashtab = rrhashtab[i];
	}
	return aux;
}*/




int Aligner::alignFree_SV(kseq_t *trunk, int offset, char *useread, Sam_Rec *_sams)
{
		
					uint32_t c_align=0, max[2];
					//int msc=0;
					uint32_t canReadLen=LEN_BASES;
					useread = useread + offset ;
					
					
				//find out the target(s)
					if(preserved.hit_times[0] >preserved.hit_times[1])
					{
						max[0]=0;
						max[1]=1;
					}
					else
					{
						max[0]=1;
						max[1]=0;
					}
				for (uint32_t ii=2;ii<canN/2 && preserved.seq_num[ii]!=0;ii++)
				{
					if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					{
						max[1]=max[0];
						max[0]=ii;
					}
					else
					{
						if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							max[1]=ii;
					}
					
				}
				
				for (uint32_t ii=canN/2;ii<canN && preserved.seq_num[ii]!=0;ii++)
				{
					if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					{
						max[1]=max[0];
						max[0]=ii;
					}
					else
					{
						if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							max[1]=ii;
					}
					
				}
				
				//cout<<"split-preserved="<<preserved.seq_num[max[0]]<<"\t"<<preserved.hit_times[max[0]]<<"\t"<<preserved.isrc[max[0]]<<endl;
				//cout<<"split-preserved="<<preserved.seq_num[max[1]]<<"\t"<<preserved.hit_times[max[1]]<<"\t"<<preserved.isrc[max[1]]<<endl;
				
				for(uint32_t ii=0;ii<2 && preserved.hit_times[max[ii]]>=1;ii++)
				{	
					
					
					if(preserved.isrc[max[ii]]==0)
						_sams[c_align].cigar="+";
					else
					    _sams[c_align].cigar="-";
					    
					_sams[c_align].score=0;
					uint8_t  ind;
					uint32_t r_startP;
					uint32_t chrstartPos;
					
					r_startP=preserved.seq_num[max[ii]];
						
						
					for (ind=1;ind<countChr;++ind)
					if (Start_pos[ind]>r_startP)
						break;
					chrstartPos = r_startP - Start_pos[ind-1];
					_sams[c_align].chrIndex = ind;
					_sams[c_align].pos = chrstartPos;
					_sams[c_align].chrLen=Start_pos[ind]-Start_pos[ind-1];
					
					
		
					
	  
			c_align++;	
			
				
			uint32_t dff=preserved.hit_times[max[0]]-preserved.hit_times[max[1]];
			if(ii==0 && dff > preserved.hit_times[max[0]]/2)
			    break;
			    
		
			}
			
			//cout<<"Split alignment\n";
		if(c_align >=1){
			
			if (c_align >= 2) {
				
				
				int quality=60;
				
				for (uint32_t i=0;i<c_align;++i) {
					if(i!=0)
					{
					
					_sams[i].MAQ = 0;
					}
					else
						_sams[i].MAQ = quality;

					
				}
				//may be if shoud be added
				
			}
		  else  {

		 _sams[0].MAQ = 60;

		}
		
		// write PAF
		
		
	 for(int i=0;i<c_align;++i)
	{
		
		
	 cout<<trunk->name.s<<"\t"<<trunk->seq.l<<"\t"<<offset<<"\t"<<offset+LEN_BASES-1<<"\t"<<_sams[i].cigar<<"\t"<<ChrName[_sams[i].chrIndex]<<"\t"<<Start_pos[_sams[i].chrIndex]-Start_pos[_sams[i].chrIndex-1]<<"\t"<<_sams[i].pos<<"\t"<<_sams[i].pos+LEN_BASES-1<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<_sams[i].MAQ<<endl;
						
						
	}
  }
		return(c_align);		
} 




int Aligner::align_SV(kseq_t *trunk, int offset, char *useread, Sam_Rec *_sams)
{
		
					uint32_t c_align=0, max[2];
					//int msc=0;
					uint32_t canReadLen=LEN_BASES;
					useread = useread + offset ;
					
					
				//find out the target(s)
					if(preserved.hit_times[0] >preserved.hit_times[1])
					{
						max[0]=0;
						max[1]=1;
					}
					else
					{
						max[0]=1;
						max[1]=0;
					}
				for (uint32_t ii=2;ii<canN;ii++)
				{
					if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					{
						max[1]=max[0];
						max[0]=ii;
					}
					else
					{
						if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							max[1]=ii;
					}
					
				}
				
				//cout<<"preserved="<<preserved.seq_num[max[0]]<<"\t"<<preserved.hit_times[max[0]]<<"\t"<<preserved.isrc[max[0]]<<endl;
				//cout<<"preserved="<<preserved.seq_num[max[1]]<<"\t"<<preserved.hit_times[max[1]]<<"\t"<<preserved.isrc[max[1]]<<endl;
				
				for(uint32_t ii=0;ii<2 && preserved.hit_times[max[ii]]>=1;ii++)
				{	
					
					
					_sams[c_align].cigar="";
					_sams[c_align].score=0;
					int n_cigar = 0;
					int qlen = 0;
					int tlen = 0;
					uint32_t *cigar;
					const 	uint8_t *readqry_;
					const 	uint8_t *refqry_;
					uint8_t  ind;
					uint32_t r_startP;
					uint32_t chrstartPos;
					const 	char 	correspondTable[] = "MIDNSHP=X";
					//uint32_t *chrStartP
					//int h0 = 0;
					uint8_t			readqry[canReadLen];
					uint8_t 		refqry[canReadLen];
					transIntoDec(readqry,useread, canReadLen);
					transIntoDec(refqry,genome+preserved.seq_num[max[ii]],canReadLen);
					readqry_ = readqry;
					refqry_ = refqry;
					
					
					//_sams[c_align].score = ksw_global2(canReadLen, readqry_, canReadLen, refqry_, 5, mat, 0, 0, 0, 0, 40, &n_cigar, &cigar);
				    _sams[c_align].score = ksw_extend_core(canReadLen, readqry_, canReadLen, refqry_, 5, mat, gapopen, gapextend, 40, canReadLen, &qlen, &tlen, &cigar, &n_cigar);
				    //cout << n_cigar<<"\n";
				    	//cout<<"here1\n";
					r_startP=preserved.seq_num[max[ii]];
						
						
					for (ind=1;ind<countChr;++ind)
					if (Start_pos[ind]>r_startP)
						break;
					chrstartPos = r_startP - Start_pos[ind-1];
					_sams[c_align].chrIndex = ind;
					_sams[c_align].pos = chrstartPos;
					
					if (preserved.isrc[max[ii]])
						_sams[c_align].flag = 16;
					else
						_sams[c_align].flag = 0;
					
							int startPosCigar = 0;
							
		char  	 trans_cigar[50];
		uint32_t  countM = 0;
		//if (qlen != canReadLen ) {
		//startPosCigar += sprintf(trans_cigar, "%uS", canReadLen - qlen);
		//_sams[c_align].cigar.append(trans_cigar);
		//}// proves that softclipings do exist

		if (n_cigar) {
		for (int z=n_cigar-1;z>0;--z) {
			
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
			_sams[c_align].cigar.append(trans_cigar);
			//++startPosCigar;
		}
		
		if (correspondTable[cigar[0]&0xf] == 'M')
			countM = cigar[0] >> 4;
		else {
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
			_sams[c_align].cigar.append(trans_cigar);
		}

		free(cigar);
		}
					
	  c_align++;
			
			
				
				uint32_t dff=preserved.hit_times[max[0]]-preserved.hit_times[max[1]];
			if(ii==0 && dff > preserved.hit_times[max[0]]/2)
			    break;
			
				
			
			}
			
			//cout<<"Split No_of alignment="<<c_align<<"\n";
		if(c_align >=1){
			
			 //cout<<"split alignment\n";
			if (c_align >= 2) {
				int best, sbest, bp;
				best=_sams[0].score;
				sbest=0;
				bp=0;
				for(int i=1;i<c_align;++i)
				{
					if(_sams[i].score>best)
					{
						sbest=best;
						best=_sams[i].score;
						bp=i;
					}
					else if(_sams[i].score>sbest)
						sbest=_sams[i].score;
					
					
				}
				//qsort_r(orders,countSam,sizeof(int),compare_sam,_sams);
				int quality;
				if (best > 0) {
					quality = (int)(250.0 * 0.25 * (double)(best- sbest)/(double)(best));
					if (quality>=60) quality = 60;
				} else quality = 0;
				
				for (uint32_t i=0;i<c_align;++i) {
					if(i!=bp)
					{
					_sams[i].flag += 256;
					_sams[i].MAQ = 0;
					}
					else
						_sams[i].MAQ = quality;

					
				}
				//may be if shoud be added
				
			}
		  else  {

		 _sams[0].MAQ = 60;

		}
		
		// write SAM
		
		char *probQual = "*";
		char *usedqual;
	 for(int i=0;i<c_align;++i)
	{
		
						if (trunk->qual.s == NULL)// this might happend so qual is wrong?
							usedqual = probQual;
						else
							usedqual = trunk->qual.s;

	 cout<<trunk->name.s<<"\t"<<_sams[i].flag<<"\t"<<ChrName[_sams[i].chrIndex]<<"\t"<<_sams[i].pos<<"\t"<<_sams[i].MAQ<<"\t"<<_sams[i].cigar<<"\t"<<"*"<<"\t"<<"0"<<"\t"
						<<"0"<<"\t";
						cout<<useread<<"\t";
						cout<<usedqual<<"\t"<<"AS:i:"<<_sams[i].score<<endl;
						
						
	}
  }
		return(c_align);		
} 


int Aligner::applySV(kseq_t *trunk, Sam_Rec *_svsams, uint32_t L, uint32_t *pat, uint32_t len_genome, int align, int tr)
{
	//split into several pieces
	uint32_t canReadLen=trunk->seq.l;
	int LEN_BS;
	
	if(canReadLen<=LEN_BASES)
		LEN_BS=LEN_BASES/2;
	else if(canReadLen<= 10*LEN_BASES)
	 LEN_BS=LEN_BASES;
	else
	 LEN_BS=2*LEN_BASES;
	
	uint32_t pNumber;
	pNumber = trunk->seq.l/LEN_BS;
	
	uint32_t leftLen =LEN_BS + trunk->seq.l - pNumber*LEN_BS;
	//if (leftLen) ++pNumber;
	// produce the
	int bkt_index[pNumber+1], max[2];
	uint32_t len[pNumber+1], amount;
	//uint32_t bulk_offset = groupNum * opt->len_limit;
	uint32_t count_aligned;
	bkt_index[0] = 0;
	len[0] = 0;
	for (uint32_t i=0;i<pNumber -1;++i) len[i+1] = len[i] + LEN_BS;
	len[pNumber] = len[pNumber-1] + leftLen;

	char *useread;
	bool isRC;
	int  movement , is_aligned=0;
	
		//Intialize preserved array
		//preserved=new bkt2;
		preserved.seq_num=new uint32_t[canN]();
		preserved.hit_times=new uint16_t[canN]();
		preserved.isrc=new bool[canN];
	
	
	if(tr==4)
	amount=((LEN_BS - 32)>>2)*L ;
	else if(tr==8)
	amount=((LEN_BS - 32)>>3)*L ;// to divide by 8 as len_sed=32
	else if(tr==16)
	 amount=((LEN_BS - 32)>>4)*L;
	 else if (tr==32)
	 amount=((LEN_BS - 32)>>5)*L;



	for (uint32_t ij = 0; ij < pNumber - 1;++ij) {    // for each split
		useread = trunk->seq.s;
		isRC = false;
		//int c;
		movement = ij*LEN_BS;
		
		
		uint32_t *rhv, chv, ii, flg, end, start;
		
		uint32_t kk, pos;
		
		
		
		//cout<<"movement="<<movement<<"\n";
		for (int j=0;j<=1;++j) {
		
		
		
		rhv= new uint32_t[amount+1];
		//compute the S-conLSH from read
		chv=generate_conLSH_read(rhv, LEN_BS, useread+movement, L,pat,tr);
		
		
							if(isRC)
							 {
								start=canN/2;
								end=canN;
							  }	
							else
							{
								start=0;
								end=canN/2;
							}
		
		for(uint32_t i=0;i<chv;i++)
		{
					
					
					if(h_index->pos[rhv[i]]!=0 && h_index->j[rhv[i]]!=0)  //successful search
					{
					  
						
						for(kk=0; kk<=h_index->j[rhv[i]] ; kk++)
						{
							 
								
							 flg=0;
							  ii=start;
							pos=h_index->pos[rhv[i]]-1;
							while(ii<end && preserved.seq_num[ii]!=0) 
							{
								//there are already preserved sequences
								uint32_t diff=hv[pos+kk]>preserved.seq_num[ii] ? (hv[pos+kk]-preserved.seq_num[ii]) : (preserved.seq_num[ii]-hv[pos+kk]);
								if(diff < LEN_BS) //localisized seq is already there
								{
									//if(hv[pos+kk]<preserved.seq_num[ii]) //store the starting position of the target
										//preserved.seq_num[ii]=hv[pos+kk];
									preserved.hit_times[ii]++;
									flg=1;
									break;
								}
								else if(diff < canReadLen)
								{
									preserved.hit_times[ii]++;
								} 
								
								ii++;
							}
							
							if(flg==0 && ii<end) // insert the new seq_num here
							{
								preserved.seq_num[ii]=hv[pos+kk];
								preserved.hit_times[ii]=1;
								preserved.isrc[ii]=isRC;
							}
							
						}
							

						//cout<<kk<<endl;
						//total_n+=s->no_of_seqs;
						
					}
				
			       // else
						
		}
			
			
			//cout<<"\n Preserved Array before=\n";
		 //for(int ll=start;ll<end && preserved.seq_num[ll]!=0;ll++)
				//cout<<preserved.seq_num[ll]<<"\t"<<preserved.hit_times[ll]<<"\t"<<preserved.isrc[ll]<<endl;
		
			
			uint32_t nm=start, kl=start;
			while(nm<end && preserved.seq_num[nm]!=0)
			{
				//kl points location of smaller hits to flush, nm points next valuable hit location
				
				while(kl<end && preserved.seq_num[kl]!=0 && preserved.hit_times[kl]>L)
					kl++;
				if(nm==start)
					nm=kl+1;
					
				while(nm<end && preserved.seq_num[nm]!=0 && preserved.hit_times[nm]<=L)
					nm++;
					
				if(nm<end && preserved.seq_num[nm]!=0)
				{
					preserved.seq_num[kl]=preserved.seq_num[nm];
					preserved.hit_times[kl]=preserved.hit_times[nm];
					preserved.isrc[kl]=preserved.isrc[nm];
					
					preserved.hit_times[nm]=0;
					kl++;
				}
				
					
			}
			
			while(kl<=nm)
			{
				preserved.seq_num[kl]=0;
				kl++;
			}
			
			
			
			//if(preserved.hit_times[start+ij*2] >preserved.hit_times[start+ij*2+1])
					//{
						//max[0]=start+ij*2;
						//max[1]=start+ij*2+1;
					//}
					//else
					//{
						//max[0]=start+ij*2+1;
						//max[1]=start+ij*2;
					//}
				//for (uint32_t ii=start+ij*2+2;ii<end && preserved.seq_num[ii]!=0;ii++)
				//{
					//if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					//{
						//max[1]=max[0];
						//max[0]=ii;
					//}
					//else
					//{
						//if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							//max[1]=ii;
						//else
						 //preserved.seq_num[ii]=0;
					//}
					
				//}
				
				////cout<<"preserved="<<preserved.seq_num[max[0]]<<"\t"<<preserved.hit_times[max[0]]<<"\t"<<preserved.isrc[max[0]]<<endl;
				////cout<<"preserved="<<preserved.seq_num[max[1]]<<"\t"<<preserved.hit_times[max[1]]<<"\t"<<preserved.isrc[max[1]]<<endl;
				////cout<<max[0]<<"\t"<<max[1]<<endl;
				
				
				//if(start+ij*2!=max[1])
				//{
				//preserved.seq_num[start+ij*2]=preserved.seq_num[max[0]];
				//preserved.hit_times[start+ij*2]=preserved.hit_times[max[0]];
				//preserved.isrc[start+ij*2]=preserved.isrc[max[0]];
				
				//preserved.seq_num[start+ij*2+1]=preserved.seq_num[max[1]];
				//preserved.hit_times[start+ij*2+1]=preserved.hit_times[max[1]];
				//preserved.isrc[start+ij*2+1]=preserved.isrc[max[1]];
				
			    //}
				//else if(start+ij*2+1!=max[0])
				//{
				//preserved.seq_num[start+ij*2+1]=preserved.seq_num[max[1]];
				//preserved.hit_times[start+ij*2+1]=preserved.hit_times[max[1]];
				//preserved.isrc[start+ij*2+1]=preserved.isrc[max[1]];
				
				//preserved.seq_num[start+ij*2]=preserved.seq_num[max[0]];
				//preserved.hit_times[start+ij*2]=preserved.hit_times[max[0]];
				//preserved.isrc[start+ij*2]=preserved.isrc[max[0]];
			    
			    //}
			    
			    
				
				//if((start+ij*2!=max[0])&&(start+ij*2+1!=max[0]))
						//preserved.seq_num[max[0]]=0;
						
				//if((start+ij*2+1!=max[1])&&(start+ij*2!=max[1]))
						//preserved.seq_num[max[1]]=0;
				
				
				
				
				
				
				
			//cout<<"\n Preserved Array After=\n";
		 //for(int ll=start;ll<end && preserved.seq_num[ll]!=0;ll++)
				//cout<<preserved.seq_num[ll]<<"\t"<<preserved.hit_times[ll]<<"\t"<<preserved.isrc[ll]<<endl;
		
			
			useread = trunk->seq.rs;
			movement = trunk->seq.l - len[ij+1];
			isRC = true;
			
			
			delete rhv;	
			
		}
		
			
		
		
		
	}


		//if(align==1)
		   //count_aligned=align_NonSV(trunk, useread,_svsams);
		//else
		   //count_aligned=alignFree_NonSV(trunk,useread,_svsams);
		   
			//if(count_aligned)
				//is_aligned=1;
		

			
			//delete preserved.seq_num;
			//delete preserved.hit_times;
			//delete preserved.isrc;


	//deal with the last piece
	
	
	//uint32_t offset = ((leftLen - LEN_BASES)>>1);
	useread = trunk->seq.s;
	isRC = false;
	movement = (pNumber-1)*LEN_BS;
	 uint32_t ij=pNumber-1;
	
		uint32_t *rhv, chv, ii, flg, end, start;
		
		uint32_t kk, pos;
		
	
	//int c;
	for (int j=0;j<=1;++j) {
		
		
		
		//amount=((LEN_BASES - 32)>>4)*L ;// to divide by 32 as len_sed=32
		
		
		
		rhv= new uint32_t[amount+1];
		//compute the S-conLSH from read
		chv=generate_conLSH_read(rhv, LEN_BS, useread+movement, L,pat,tr);
		
		
							if(isRC)
							 {
								start=canN/2;
								end=canN;
							  }	
							else
							{
								start=0;
								end=canN/2;
							}
		
		for(uint32_t i=0;i<chv;i++)
		{
					
					
					if(h_index->pos[rhv[i]]!=0 && h_index->j[rhv[i]]!=0)  //successful search
					{
					  
						
						for(kk=0; kk<=h_index->j[rhv[i]] ; kk++)
						{
							 
								
							 flg=0;
							  ii=start;
							pos=h_index->pos[rhv[i]]-1;
							while(ii<end && preserved.seq_num[ii]!=0) 
							{
								//there are already preserved sequences
								uint32_t diff=hv[pos+kk]>preserved.seq_num[ii] ? (hv[pos+kk]-preserved.seq_num[ii]) : (preserved.seq_num[ii]-hv[pos+kk]);
								if(diff < LEN_BS) //localisized seq is already there
								{
									//if(hv[pos+kk]<preserved.seq_num[ii]) //store the starting position of the target
										//preserved.seq_num[ii]=hv[pos+kk];
									preserved.hit_times[ii]++;
									flg=1;
									break;
								}
								else if(diff < canReadLen)
								{
									preserved.hit_times[ii]++;
								} 
								
								ii++;
							}
							
							if(flg==0 && ii<end) // insert the new seq_num here
							{
								preserved.seq_num[ii]=hv[pos+kk];
								preserved.hit_times[ii]=1;
								preserved.isrc[ii]=isRC;
							}
							
						}
							

						//cout<<kk<<endl;
						//total_n+=s->no_of_seqs;
						
					}
				
			       // else
						
		}
		
		
		
				//if(preserved.hit_times[start+ij*2] >preserved.hit_times[start+ij*2+1])
					//{
						//max[0]=start+ij*2;
						//max[1]=start+ij*2+1;
					//}
					//else
					//{
						//max[0]=start+ij*2+1;
						//max[1]=start+ij*2;
					//}
				//for (uint32_t ii=start+ij*2+2;ii<end && preserved.seq_num[ii]!=0;ii++)
				//{
					//if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					//{
						//max[1]=max[0];
						//max[0]=ii;
					//}
					//else
					//{
						//if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							//max[1]=ii;
						//else
						 //preserved.seq_num[ii]=0;
					//}
					
				//}
				
			
			 
				//if(start+ij*2!=max[1])
				//{
				//preserved.seq_num[start+ij*2]=preserved.seq_num[max[0]];
				//preserved.hit_times[start+ij*2]=preserved.hit_times[max[0]];
				//preserved.isrc[start+ij*2]=preserved.isrc[max[0]];
				
				//preserved.seq_num[start+ij*2+1]=preserved.seq_num[max[1]];
				//preserved.hit_times[start+ij*2+1]=preserved.hit_times[max[1]];
				//preserved.isrc[start+ij*2+1]=preserved.isrc[max[1]];
				
			    //}
				//else if(start+ij*2+1!=max[0])
				//{
				//preserved.seq_num[start+ij*2+1]=preserved.seq_num[max[1]];
				//preserved.hit_times[start+ij*2+1]=preserved.hit_times[max[1]];
				//preserved.isrc[start+ij*2+1]=preserved.isrc[max[1]];
				
				//preserved.seq_num[start+ij*2]=preserved.seq_num[max[0]];
				//preserved.hit_times[start+ij*2]=preserved.hit_times[max[0]];
				//preserved.isrc[start+ij*2]=preserved.isrc[max[0]];
			    
			    //}
			    
			    
				
				//if((start+ij*2!=max[0])&&(start+ij*2+1!=max[0]))
						//preserved.seq_num[max[0]]=0;
						
				//if((start+ij*2+1!=max[1])&&(start+ij*2!=max[1]))
						//preserved.seq_num[max[1]]=0;
				
				
			//cout<<"\n Preserved Array=\n";
		 //for(int ll=start;ll<end && preserved.seq_num[ll]!=0;ll++)
				//cout<<preserved.seq_num[ll]<<"\t"<<preserved.hit_times[ll]<<"\t"<<preserved.isrc[ll]<<endl;
		
		
		
		//cout<<"\n Preserved Array before=\n";
		 //for(int ll=start;ll<end && preserved.seq_num[ll]!=0;ll++)
				//cout<<preserved.seq_num[ll]<<"\t"<<preserved.hit_times[ll]<<"\t"<<preserved.isrc[ll]<<endl;
		
			
			uint32_t nm=start, kl=start;
			while(nm<end && preserved.seq_num[nm]!=0)
			{
				//kl points location of smaller hits to flush, nm points next valuable hit location
				
				while(kl<end && preserved.seq_num[kl]!=0 && preserved.hit_times[kl]>L)
					kl++;
				if(nm==start)
					nm=kl+1;
					
				while(nm<end && preserved.seq_num[nm]!=0 && preserved.hit_times[nm]<=L)
					nm++;
					
				if(nm<end && preserved.seq_num[nm]!=0)
				{
					preserved.seq_num[kl]=preserved.seq_num[nm];
					preserved.hit_times[kl]=preserved.hit_times[nm];
					preserved.isrc[kl]=preserved.isrc[nm];
					
					preserved.hit_times[nm]=0;
					kl++;
				}
				
					
			}
			
			while(kl<=nm)
			{
				preserved.seq_num[kl]=0;
				kl++;
			}
			
			
			
			//cout<<"\n Preserved Array before=\n";
		 //for(int ll=start;ll<end && preserved.seq_num[ll]!=0;ll++)
				//cout<<preserved.seq_num[ll]<<"\t"<<preserved.hit_times[ll]<<"\t"<<preserved.isrc[ll]<<endl;
		
		
		useread = trunk->seq.rs ;
		movement = 0;
		isRC = true;
		
		delete rhv;
	}

		////if(align==1)
		   //count_aligned=align_SV(trunk,movement, useread,_svsams);
		//else
		   //count_aligned=alignFree_SV(trunk,movement, useread,_svsams);
		
			 //if(count_aligned)
				//is_aligned=1;
		
			//cout<<"\n Preserved Array appplySV=\n";
		// for(int ll=start;ll<end && preserved.seq_num[ll]!=0;ll++)
				//cout<<preserved.seq_num[ll]<<"\t"<<preserved.hit_times[ll]<<"\t"<<preserved.isrc[ll]<<endl;
		if(align==0)
		 count_aligned=alignFree_NonSV(trunk,useread,_svsams);
		 else
		 count_aligned=align_NonSV(trunk,useread,_svsams);
		 
		 	
			delete preserved.seq_num;
			delete preserved.hit_times;
			delete preserved.isrc;


  return(is_aligned);
}


//int Aligner::alignFree_NonSV(kseq_t *trunk, char *useread ,Sam_Rec *_sams)
//{
					//uint32_t c_align=0, max[2], chim=0;
					////int msc=0;
					//uint32_t canReadLen=trunk->seq.l;
					////char * useread = trunk->seq.s;
					
					////find out the target(s)
					//if(preserved.hit_times[0] >preserved.hit_times[1])
					//{
						//max[0]=0;
						//max[1]=1;
					//}
					//else
					//{
						//max[0]=1;
						//max[1]=0;
					//}
				//for (uint32_t ii=2;ii<canN/2 && preserved.seq_num[ii]!=0;ii++)
				//{
					//if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					//{
						//max[1]=max[0];
						//max[0]=ii;
					//}
					//else
					//{
						//if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							//max[1]=ii;
					//}
					
				//}
				
				//for (uint32_t ii=canN/2;ii<canN && preserved.seq_num[ii]!=0;ii++)
				//{
					//if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					//{
						//max[1]=max[0];
						//max[0]=ii;
					//}
					//else
					//{
						//if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							//max[1]=ii;
					//}
					
				//}
				
				////cout<<"preserved="<<preserved.seq_num[max[0]]<<"\t"<<preserved.hit_times[max[0]]<<"\t"<<preserved.isrc[max[0]]<<endl;
				////cout<<"preserved="<<preserved.seq_num[max[1]]<<"\t"<<preserved.hit_times[max[1]]<<"\t"<<preserved.isrc[max[1]]<<endl;
				
				//for(uint32_t ii=0;ii<2 &&  preserved.hit_times[max[ii]]>=1;ii++)
				//{
					//uint32_t df=preserved.hit_times[max[0]]-preserved.hit_times[max[1]];
					//if(((df < preserved.hit_times[max[0]]/2) || (df ==0)) && (canReadLen > LEN_BASES) )
						//chim=1;
				////&& ((preserved.hit_times[max[0]]<10) || ((preserved.hit_times[max[0]]/(float)canReadLen)<0.003) )
					//if(preserved.isrc[max[ii]]==0)
						//_sams[c_align].cigar="+";
					//else
					    //_sams[c_align].cigar="-";
					    
					//_sams[c_align].score=0;
					//uint8_t  ind;
					//uint32_t r_startP;
					//uint32_t chrstartPos;
					
					//r_startP=preserved.seq_num[max[ii]];
						
						
					//for (ind=1;ind<countChr;++ind)
					//if (Start_pos[ind]>r_startP)
						//break;
					//chrstartPos = r_startP - Start_pos[ind-1];
					//_sams[c_align].chrIndex = ind;
					//_sams[c_align].pos = chrstartPos;
					//_sams[c_align].chrLen=Start_pos[ind]-Start_pos[ind-1];
					
					
		
					
	  
			//c_align++;	
				
			//uint32_t dff=preserved.hit_times[max[0]]-preserved.hit_times[max[1]];
			//if(ii==0 && dff > preserved.hit_times[max[0]]/2)
			    //break;
			    
		
				
     //}
			
		//if(c_align >=1){
			
			 
			//if (c_align >= 2) {
				
				
				//int quality=60;
				
				//for (uint32_t i=0;i<c_align;++i) {
					//if(i!=0)
					//{
					
					//_sams[i].MAQ = 0;
					//}
					//else
						//_sams[i].MAQ = quality;

					
				//}
				////may be if shoud be added
				
			//}
		  //else  {

		 //_sams[0].MAQ = 60;

		//}
		
		//// write PAF
		
		
	 //for(int i=0;i<c_align;++i)
	//{
		
		
	 //cout<<trunk->name.s<<"\t"<<canReadLen<<"\t"<<"0"<<"\t"<<canReadLen-1<<"\t"<<_sams[i].cigar<<"\t"<<ChrName[_sams[i].chrIndex]<<"\t"<<Start_pos[_sams[i].chrIndex]-Start_pos[_sams[i].chrIndex-1]<<"\t"<<_sams[i].pos<<"\t"<<_sams[i].pos+canReadLen-1<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<_sams[i].MAQ<<endl;
						
						
	//}
  //} 
      ////delete _sams;
      //if(chim==1)
		//return(0);
	 //else
		//return(c_align);		
//}






int Aligner::alignFree_NonSV(kseq_t *trunk, char *useread ,Sam_Rec *_sams)
{
					uint32_t c_align=0, max[2], chim=0;
					//int msc=0;
					uint32_t canReadLen=trunk->seq.l;
					//char * useread = trunk->seq.s;
					
					//find out the target(s)
					if(preserved.hit_times[0] >preserved.hit_times[1])
					{
						max[0]=0;
						max[1]=1;
					}
					else
					{
						max[0]=1;
						max[1]=0;
					}
				for (uint32_t ii=2;ii<canN/2 && preserved.seq_num[ii]!=0;ii++)
				{
					if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					{
						max[1]=max[0];
						max[0]=ii;
					}
					else
					{
						if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							max[1]=ii;
					}
					
				}
				
				for (uint32_t ii=canN/2;ii<canN && preserved.seq_num[ii]!=0;ii++) 
				{
					if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					{
						max[1]=max[0];
						max[0]=ii;
					}
					else
					{
						if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							max[1]=ii;
					}
					
				}
				
				//cout<<"preserved="<<preserved.seq_num[max[0]]<<"\t"<<preserved.hit_times[max[0]]<<"\t"<<preserved.isrc[max[0]]<<endl;
				//cout<<"preserved="<<preserved.seq_num[max[1]]<<"\t"<<preserved.hit_times[max[1]]<<"\t"<<preserved.isrc[max[1]]<<endl;
				
				for(uint32_t ii=0;ii<canN/2 && preserved.seq_num[ii]!=0;ii++)
				{
					
					
					
					uint32_t df=preserved.hit_times[max[0]]-preserved.hit_times[ii];
					//if(((df < preserved.hit_times[max[0]]/2) || (df ==0)) && (canReadLen > LEN_BASES) )
						//chim=1;
				////&& ((preserved.hit_times[max[0]]<10) || ((preserved.hit_times[max[0]]/(float)canReadLen)<0.003) )
				
				if(df < preserved.hit_times[max[0]]/2)
				{
					if(preserved.isrc[ii]==0)
						_sams[c_align].cigar="+";
					else
					    _sams[c_align].cigar="-";
					    
					_sams[c_align].score=0;
					uint8_t  ind;
					uint32_t r_startP;
					uint32_t chrstartPos;
					
					r_startP=preserved.seq_num[ii];
						
						
					for (ind=1;ind<countChr;++ind)
					if (Start_pos[ind]>r_startP)
						break;
					chrstartPos = r_startP - Start_pos[ind-1];
					_sams[c_align].chrIndex = ind;
					_sams[c_align].pos = chrstartPos;
					_sams[c_align].chrLen=Start_pos[ind]-Start_pos[ind-1];
					
					
		
					
	  
			c_align++;	
				
			//uint32_t dff=preserved.hit_times[max[0]]-preserved.hit_times[max[1]];
			//if(ii==0 && dff > preserved.hit_times[max[0]]/2)
			    //break;
		}	    
		
	}	
		for(uint32_t ii=canN/2;ii<canN && preserved.seq_num[ii]!=0;ii++)
				{
					
					
					
					uint32_t df=preserved.hit_times[max[0]]-preserved.hit_times[ii];
					//if(((df < preserved.hit_times[max[0]]/2) || (df ==0)) && (canReadLen > LEN_BASES) )
						//chim=1;
				////&& ((preserved.hit_times[max[0]]<10) || ((preserved.hit_times[max[0]]/(float)canReadLen)<0.003) )
				
				if(df < preserved.hit_times[max[0]]/2)
				{
					if(preserved.isrc[ii]==0)
						_sams[c_align].cigar="+";
					else
					    _sams[c_align].cigar="-";
					    
					_sams[c_align].score=0;
					uint8_t  ind;
					uint32_t r_startP;
					uint32_t chrstartPos;
					
					r_startP=preserved.seq_num[ii];
						
						
					for (ind=1;ind<countChr;++ind)
					if (Start_pos[ind]>r_startP)
						break;
					chrstartPos = r_startP - Start_pos[ind-1];
					_sams[c_align].chrIndex = ind;
					_sams[c_align].pos = chrstartPos;
					_sams[c_align].chrLen=Start_pos[ind]-Start_pos[ind-1];
					
					
		
					
	  
			c_align++;	
				
			//uint32_t dff=preserved.hit_times[max[0]]-preserved.hit_times[max[1]];
			//if(ii==0 && dff > preserved.hit_times[max[0]]/2)
			    //break;
		}	    
		
				
     }
			
		if(c_align >=1){
			
			 
			if (c_align >= 2) {
				
				
				int quality=60;
				
				for (uint32_t i=0;i<c_align;++i) {
					if(i!=0)
					{
					
					_sams[i].MAQ = 0;
					}
					else
						_sams[i].MAQ = quality;

					
				}
				//may be if shoud be added
				
			}
		  else  {

		 _sams[0].MAQ = 60;

		}
		
		// write PAF
		
		
	 for(int i=0;i<c_align;++i)
	{
		
		
	 cout<<trunk->name.s<<"\t"<<canReadLen<<"\t"<<"0"<<"\t"<<canReadLen-1<<"\t"<<_sams[i].cigar<<"\t"<<ChrName[_sams[i].chrIndex]<<"\t"<<Start_pos[_sams[i].chrIndex]-Start_pos[_sams[i].chrIndex-1]<<"\t"<<_sams[i].pos<<"\t"<<_sams[i].pos+canReadLen-1<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<_sams[i].MAQ<<endl;
						
						
	}
  } 
      //delete _sams;
      if(chim==1)
		return(0);
	 else
		return(c_align);		
}








int Aligner::align_NonSV(kseq_t *trunk, char *useread ,Sam_Rec *_sams)
{
				
			uint32_t c_align=0, max[2], chim=0;
					//int msc=0;
					uint32_t canReadLen=trunk->seq.l;
					//char * useread = trunk->seq.s;
					
					//find out the target(s)
					if(preserved.hit_times[0] >preserved.hit_times[1])
					{
						max[0]=0;
						max[1]=1;
					}
					else
					{
						max[0]=1;
						max[1]=0;
					}
				for (uint32_t ii=2;ii<canN/2 && preserved.seq_num[ii]!=0;ii++)
				{
					if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					{
						max[1]=max[0];
						max[0]=ii;
					}
					else
					{
						if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							max[1]=ii;
					}
					
				}
				
				for (uint32_t ii=canN/2;ii<canN && preserved.seq_num[ii]!=0;ii++) 
				{
					if(preserved.hit_times[ii] >preserved.hit_times[max[0]])
					{
						max[1]=max[0];
						max[0]=ii;
					}
					else
					{
						if(preserved.hit_times[ii] > preserved.hit_times[max[1]])
							max[1]=ii;
					}
					
				}
				
				//cout<<"preserved="<<preserved.seq_num[max[0]]<<"\t"<<preserved.hit_times[max[0]]<<"\t"<<preserved.isrc[max[0]]<<endl;
				//cout<<"preserved="<<preserved.seq_num[max[1]]<<"\t"<<preserved.hit_times[max[1]]<<"\t"<<preserved.isrc[max[1]]<<endl;
				
				for(uint32_t ii=0;ii<canN/2 && preserved.seq_num[ii]!=0;ii++)
				{
					
					
					
					uint32_t df=preserved.hit_times[max[0]]-preserved.hit_times[ii];
					//if(((df < preserved.hit_times[max[0]]/2) || (df ==0)) && (canReadLen > LEN_BASES) )
						//chim=1;
				////&& ((preserved.hit_times[max[0]]<10) || ((preserved.hit_times[max[0]]/(float)canReadLen)<0.003) )
				
				if(df < preserved.hit_times[max[0]]/2)
				{
					
					_sams[c_align].cigar="";
					_sams[c_align].score=0;
					int n_cigar = 0;
					int qlen = 0;
					int tlen = 0;
					uint32_t *cigar;
					const 	uint8_t *readqry_;
					const 	uint8_t *refqry_;
					uint8_t  ind;
					uint32_t r_startP;
					uint32_t chrstartPos;
					const 	char 	correspondTable[] = "MIDNSHP=X";
					
					uint8_t			readqry[canReadLen];
					uint8_t 		refqry[canReadLen];
					transIntoDec(readqry,useread, canReadLen);
					transIntoDec(refqry,genome+preserved.seq_num[ii],canReadLen);
					readqry_ = readqry;
					refqry_ = refqry;
					
					
					//_sams[c_align].score = ksw_global2(canReadLen, readqry_, canReadLen, refqry_, 5, mat, 0, 0, 0, 0, 40, &n_cigar, &cigar);
				    _sams[c_align].score = ksw_extend_core(canReadLen, readqry_, canReadLen, refqry_, 5, mat, gapopen, gapextend, 40, canReadLen, &qlen, &tlen, &cigar, &n_cigar);
				    //cout << n_cigar<<"\n";
				    	//cout<<"here1\n";
					r_startP=preserved.seq_num[ii];
						
						
					for (ind=1;ind<countChr;++ind)
					if (Start_pos[ind]>r_startP)
						break;
					chrstartPos = r_startP - Start_pos[ind-1];
					_sams[c_align].chrIndex = ind;
					_sams[c_align].pos = chrstartPos;
					
					if (preserved.isrc[ii])
						_sams[c_align].flag = 16;
					else
						_sams[c_align].flag = 0;
					
		int startPosCigar = 0;
		char  	 trans_cigar[50];
		uint32_t  countM = 0;
		//if (qlen != canReadLen ) {
		//startPosCigar += sprintf(trans_cigar, "%uS", canReadLen - qlen);
		//_sams[c_align].cigar.append(trans_cigar);
		//}// proves that softclipings do exist

		if (n_cigar) {
		for (int z=n_cigar-1;z>0;--z) {
			
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
			_sams[c_align].cigar.append(trans_cigar);
			//++startPosCigar;
		}
		
		if (correspondTable[cigar[0]&0xf] == 'M')
			countM = cigar[0] >> 4;
		else {
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
			_sams[c_align].cigar.append(trans_cigar);
		}

		free(cigar);
		}
					
	  c_align++;
				
				
			//uint32_t dff=preserved.hit_times[max[0]]-preserved.hit_times[max[1]];
			//if(ii==0 && dff > preserved.hit_times[max[0]]/2)
			    //break;
		}	    
		
	}	
		for(uint32_t ii=canN/2;ii<canN && preserved.seq_num[ii]!=0;ii++)
				{
					
					
					
					uint32_t df=preserved.hit_times[max[0]]-preserved.hit_times[ii];
					//if(((df < preserved.hit_times[max[0]]/2) || (df ==0)) && (canReadLen > LEN_BASES) )
						//chim=1;
				////&& ((preserved.hit_times[max[0]]<10) || ((preserved.hit_times[max[0]]/(float)canReadLen)<0.003) )
				
				if(df < preserved.hit_times[max[0]]/2)
				{
					_sams[c_align].cigar="";
					_sams[c_align].score=0;
					int n_cigar = 0;
					int qlen = 0;
					int tlen = 0;
					uint32_t *cigar;
					const 	uint8_t *readqry_;
					const 	uint8_t *refqry_;
					uint8_t  ind;
					uint32_t r_startP;
					uint32_t chrstartPos;
					const 	char 	correspondTable[] = "MIDNSHP=X";
					//uint32_t *chrStartP
					//int h0 = 0;
					uint8_t			readqry[canReadLen];
					uint8_t 		refqry[canReadLen];
					transIntoDec(readqry,useread, canReadLen);
					transIntoDec(refqry,genome+preserved.seq_num[ii],canReadLen);
					readqry_ = readqry;
					refqry_ = refqry;
					
					
					//_sams[c_align].score = ksw_global2(canReadLen, readqry_, canReadLen, refqry_, 5, mat, 0, 0, 0, 0, 40, &n_cigar, &cigar);
				    _sams[c_align].score = ksw_extend_core(canReadLen, readqry_, canReadLen, refqry_, 5, mat, gapopen, gapextend, 40, canReadLen, &qlen, &tlen, &cigar, &n_cigar);
				    //cout << n_cigar<<"\n";
				    	//cout<<"here1\n";
					r_startP=preserved.seq_num[ii];
						
						
					for (ind=1;ind<countChr;++ind)
					if (Start_pos[ind]>r_startP)
						break;
					chrstartPos = r_startP - Start_pos[ind-1];
					_sams[c_align].chrIndex = ind;
					_sams[c_align].pos = chrstartPos;
					
					if (preserved.isrc[ii])
						_sams[c_align].flag = 16;
					else
						_sams[c_align].flag = 0;
					
		int startPosCigar = 0;
		char  	 trans_cigar[50];
		uint32_t  countM = 0;
		//if (qlen != canReadLen ) {
		//startPosCigar += sprintf(trans_cigar, "%uS", canReadLen - qlen);
		//_sams[c_align].cigar.append(trans_cigar);
		//}// proves that softclipings do exist

		if (n_cigar) {
		for (int z=n_cigar-1;z>0;--z) {
			
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
			_sams[c_align].cigar.append(trans_cigar);
			//++startPosCigar;
		}
		
		if (correspondTable[cigar[0]&0xf] == 'M')
			countM = cigar[0] >> 4;
		else {
			startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
			_sams[c_align].cigar.append(trans_cigar);
		}

		free(cigar);
		}
					
	  c_align++;
				
				
				
			//uint32_t dff=preserved.hit_times[max[0]]-preserved.hit_times[max[1]];
			//if(ii==0 && dff > preserved.hit_times[max[0]]/2)
			    //break;
		}	    
		
				
     }
		if(c_align >=1){
			
			 //cout<<"gotit\n";
			if (c_align >= 2) {
				int best, sbest, bp;
				best=_sams[0].score;
				sbest=0;
				bp=0;
				for(uint32_t i=1;i<c_align;++i)
				{
					if(_sams[i].score>best)
					{
						sbest=best;
						best=_sams[i].score;
						bp=i;
					}
					else if(_sams[i].score>sbest)
						sbest=_sams[i].score;
					
					
				}
				//qsort_r(orders,countSam,sizeof(int),compare_sam,_sams);
				int quality;
				if (best > 0) {
					quality = (int)(250.0 * 0.25 * (double)(best- sbest)/(double)(best));
					if (quality>=60) quality = 60;
				} else quality = 0;
				
				for (uint32_t i=0;i<c_align;++i) {
					if(i!=bp)
					{
					_sams[i].flag += 256;
					_sams[i].MAQ = 0;
					}
					else
						_sams[i].MAQ = quality;

					
				}
				//may be if shoud be added
				
			}
		  else  {

		 _sams[0].MAQ = 60;

		}
		
		// write SAM
		
		char *probQual = "*";
		char *usedqual;
	 for(int i=0;i<c_align;++i)
	{
		
						if (trunk->qual.s == NULL)// this might happend so qual is wrong?
							usedqual = probQual;
						else
							usedqual = trunk->qual.s;

	 cout<<trunk->name.s<<"\t"<<_sams[i].flag<<"\t"<<ChrName[_sams[i].chrIndex]<<"\t"<<_sams[i].pos<<"\t"<<_sams[i].MAQ<<"\t"<<_sams[i].cigar<<"\t"<<"*"<<"\t"<<"0"<<"\t"
						<<"0"<<"\t";
						cout<<useread<<"\t";
						cout<<usedqual<<"\t"<<"AS:i:"<<_sams[i].score<<endl;
						
						
	}
  } 
      //delete _sams;
     
		return(c_align);							
								
				
}

int  Aligner::applyNonSV(kseq_t *trunk, Sam_Rec *_sams, uint32_t L, uint32_t *pat, uint32_t len_genome, int align, int tr)
{
	

	//uint32_t diff = trunk->seq.l > LEN_BASES?trunk->seq.l - LEN_BASES:0;

	//uint32_t offset = diff>>1;

	//uint32_t canReadLen = trunk->seq.l > LEN_BASES ? LEN_BASES:trunk->seq.l;
	//cout<<trunk->seq.l<<"\t"<<canReadLen;
	uint32_t canReadLen=trunk->seq.l;
	int is_aligned=0;
	char * useread = trunk->seq.s;
	bool isRC = false;
	
	
	uint32_t *rhv, chv, ii, flg, amount, end, start;
		
		uint32_t count_aligned,kk, pos;
		
		
		//Intialize preserved array
		//preserved=new bkt2;
		preserved.seq_num=new uint32_t[canN]();
		preserved.hit_times=new uint16_t[canN]();
		preserved.isrc=new bool[canN];
		
		
		if(tr==4)
		 amount=((canReadLen - 32)>>2)*L ;
		else if(tr==8)
			amount=((canReadLen - 32)>>3)*L ;
		else if(tr==16)
		 amount=((canReadLen - 32)>>4)*L ;
		 else if(tr==32)
		 amount=((canReadLen - 32)>>5)*L ;
		 else if(tr==64)
		 amount=((canReadLen - 32)>>6)*L ; 
		 else if(tr==128)
		 amount=((canReadLen - 32)>>7)*L ; 
	
	for (int j=0;j<=1;j++) {
		
		 
		 
		rhv= new uint32_t[amount+1];
		//compute the S-conLSH from read
		chv=generate_conLSH_read(rhv, canReadLen, useread, L,pat,tr);
		
		
							if(isRC)
							 {
								start=canN/2;
								end=canN;
							  }	
							else
							{
								start=0;
								end=canN/2;
							}
		
		for(uint32_t i=0;i<chv;i++)
		{
					
					
					if(h_index->pos[rhv[i]]!=0 && h_index->j[rhv[i]]!=0)  //successful search
					{
					  
						
						for(kk=0; kk<=h_index->j[rhv[i]] ; kk++)
						{
							 
								
							 flg=0;
							  ii=start;
							pos=h_index->pos[rhv[i]]-1;
							while(ii<end && preserved.seq_num[ii]!=0) 
							{
								//there are already preserved sequences
								uint32_t diff=hv[pos+kk]>preserved.seq_num[ii] ? (hv[pos+kk]-preserved.seq_num[ii]) : (preserved.seq_num[ii]-hv[pos+kk]);
								if(diff < LEN_BASES) //localisized seq is already there
								{
									//if(hv[pos+kk]<preserved.seq_num[ii]) //store the starting position of the target
										//preserved.seq_num[ii]=hv[pos+kk];
									preserved.hit_times[ii]++;
									flg=1;
									break;
								}
								else if(diff < canReadLen)
								{
									preserved.hit_times[ii]++;
								} 
								ii++;
							}
							
							if(flg==0 && ii<end) // insert the new seq_num here
							{
								preserved.seq_num[ii]=hv[pos+kk];
								preserved.hit_times[ii]=1;
								preserved.isrc[ii]=isRC;
							}
							
						}
							

						//cout<<kk<<endl;
						//total_n+=s->no_of_seqs;
						
					}
				
			       // else
						
		} 
	 
		
		useread = trunk->seq.rs;
		isRC = true;
		delete rhv;
					//cout<<"\n";
	}
	
		//cout<<"\n Preserved Array ApplyNonSV=\n";
		// for(ii=0;ii<canN && preserved.seq_num[ii]!=0;ii++)
				//cout<<preserved.seq_num[ii]<<"\t"<<preserved.hit_times[ii]<<"\t"<<preserved.isrc[ii]<<endl;
		
		
		//align using KSW 
		if(align==1)
			count_aligned=align_NonSV(trunk, useread, _sams);
		else
			count_aligned=alignFree_NonSV(trunk, useread, _sams);	
			
			if(count_aligned)
				is_aligned=1;
		
			
			delete preserved.seq_num;
			delete preserved.hit_times;
			delete preserved.isrc;

	
	
	return(is_aligned);
	
}






void Aligner::Runtask()
{
	
	
	
	//read datas
	/*len_genome = 0;
	//read_file rd;

	ChrName = new char *[LEN];
	for (int i=0;i<100;++i) { ChrName[i] = new char[LEN];}
	//int 	 countChr;
	genome = rd.read_ref(opt->refpath,&len_genome,Start_pos,ChrName,&countChr);

	if (NULL == genome) {

		fprintf(stderr,"Fail to load Genome reference, now exit");
		exit(1);
	}
	genome_e = genome + len_genome - 1;*/
	
	
	
	
  if(algn==1)
  {
	
   	//output header
 	cout<<"@HD\tVN:"<<PACKAGE_VERSION<<endl;
	for (int i=1;i<countChr;++i) {cout<<"@SQ\tSN:"<<ChrName[i]<<"\tLN:"<<Start_pos[i]-Start_pos[i-1]<<endl;}
	cout<<"@PG\tID:"<<PACKAGE_NAME<<"\tVN:"<<PACKAGE_VERSION<<"\tCL:";
	for (int i=0;i<argc;++i) {cout<<argv[i]<<" ";}
	cout<<endl;
	
 
  }
		
	
		
	gzFile fp;
	kseq_t *trunk;
	fp = gzopen(readpath, "r");

	if (thread <= 1) {

		kseq_t *seqs = kseq_init(fp);

		Sam_Rec *sams = new Sam_Rec[canN];
		
		
		while (kseq_read(seqs)>=0) {

			revComRead(seqs->seq.s, seqs->seq.rs, seqs->seq.l);
				uint16_t sam_details;
			//cout<< seqs->name.s <<"\n";
			
			if(seqs->seq.l <= LEN_BASES)
			{
				sam_details = applyNonSV(seqs, sams, L, pat, len_genome,algn, 16);
				sam_details = applyNonSV(seqs, sams, L, pat, len_genome,algn, 8);
				sam_details = applyNonSV(seqs, sams, L, pat, len_genome,algn, 4);
				

			}
			else if(seqs->seq.l <= 4*LEN_BASES)
				sam_details = applyNonSV(seqs, sams, L, pat, len_genome,algn, 16);
			else if(seqs->seq.l <= 6*LEN_BASES)
				sam_details = applyNonSV(seqs, sams, L, pat, len_genome,algn, 32);
			else if(seqs->seq.l >= 10*LEN_BASES && seqs->seq.l <= 15*LEN_BASES)
				sam_details = applyNonSV(seqs, sams, L, pat, len_genome,algn, 64);
				
				
		if ( !(seqs->seq.l <= LEN_BASES/2)) {
				
				if(seqs->seq.l <= 4*LEN_BASES)
					sam_details = applySV(seqs, sams, L, pat, len_genome,algn,4);
				else if(seqs->seq.l <= 6*LEN_BASES)
					sam_details = applySV(seqs, sams, L, pat, len_genome,algn,8);
				else if(seqs->seq.l <= 10*LEN_BASES)
					sam_details = applySV(seqs, sams, L, pat, len_genome,algn,16);
				else
					sam_details = applySV(seqs, sams, L, pat, len_genome,algn,32);
			} 
		}
		kseq_destroy(seqs);
		
	} 
	
	
   
	if (NULL != genome)
		delete[] genome;
   
	
	
	gzclose(fp);

  
}
