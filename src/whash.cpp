#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
using namespace std;
#include "whash.h"
#include "conLSH.h"
#include  "aligner.h"


//uint8_t trans[128]= {
		//4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		//4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		//4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		//4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		//4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		//4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
		//4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		//4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4
	//};
	
	
                         
	



int Hash::write_hashfile(char *path,char *genome,uint32_t len_genome,uint32_t len_sed,uint32_t K,uint32_t L,uint32_t lambda,uint32_t zero)
{
	
	//char filepath[256] = {0};

	uint32_t ccv;
	ccv=(1<<2*K*lambda)+1;
	
	//Generate the pattern
	
	pat= new uint32_t[L];
	generate_pat(pat,K,L,lambda,zero);
		 
		 
		  hashtab = new Hashtab;
		if( NULL == hashtab) {
		fprintf(stderr,"Failed when applying for new space");
		exit(1);
	  }
	  
	    h_index = new H_index; 
	  h_index->pos= new uint32_t[ccv];
	  h_index->j= new uint8_t[ccv]();
	  if( NULL == h_index->pos) {
		fprintf(stderr,"Failed when applying for hash_index");
		exit(1);
	  }
	  if( NULL == h_index->j) {
		fprintf(stderr,"Failed when applying for hash_index");
		exit(1);
	  }
	  
	 uint32_t amount;
	 amount=(len_genome - len_sed) ;
	//hashtab->seq_no = new uint32_t[amount+1];
	
    hashtab->hash_val = new uint32_t[amount+1]();
   // hashtab->sqn = new uint32_t[amount+1];
    if( NULL == hashtab->hash_val) {
		fprintf(stderr,"Failed when applying for new space hash_val");
		exit(1);
	  }
	 /* if( NULL == hashtab->sqn) {
		fprintf(stderr,"Failed when applying for new space hash_val");
		exit(1);
	  }*/
	  
	  
    uint32_t weight;
    weight=K*lambda;
    // Generate conLSH values
    count=generate_conLSH(h_index->pos, h_index->j,hashtab->hash_val, len_genome, len_sed, genome,L, pat);
    
	
	uint32_t ssn,last,next=1, diff;
	
	
	//Second pass to store the hash values
	  hv=new uint32_t[count+1]();    // hv is bounded by the no_of_valid_sequences
	  
	  
  for (uint32_t ii=0;ii<amount;ii++)
 {
		
		
		//t.insert(hashtab->hash_val[ii],ii%(len_genome-opt->len_sed));
		//cout<<"total distinct hash values="<<t.n_vals<<"\n";
		
		
	if(hashtab->hash_val[ii]!=0)
	{	
		ssn=ii;
		//cout<<hashtab->hash_val[ii]<<"\t"<<ii<<"\t"<<h_index->pos[hashtab->hash_val[ii]]<<"\t"<<next<<"\n";
		if(h_index->pos[hashtab->hash_val[ii]]==0)  //first entry of this hash value
	  {
		if(h_index->j[hashtab->hash_val[ii]]!=0 &&  h_index->j[hashtab->hash_val[ii]]!=201)  //valid entry
		{
				  h_index->pos[hashtab->hash_val[ii]]=next; //starts from 1, subtract 1 at the time of new insertion
				  next+=h_index->j[hashtab->hash_val[ii]];
				  
				 // if(h_index->pos[hashtab->hash_val[ii]]>count+1)
				       // cout<<"overflow\n";
			
			
			h_index->j[hashtab->hash_val[ii]]=0;
			hv[h_index->pos[hashtab->hash_val[ii]]+h_index->j[hashtab->hash_val[ii]]-1]=ssn ;  
	     }
			
	   }
	  else
	  { 
		last=hv[h_index->pos[hashtab->hash_val[ii]]+h_index->j[hashtab->hash_val[ii]]-1];
		diff=last>ssn ? (last-ssn) : (ssn-last);
		
			if(diff < LEN_BASES)
			{
				if(last>ssn)
				{
					hv[h_index->pos[hashtab->hash_val[ii]]+h_index->j[hashtab->hash_val[ii]]-1]=ssn;
					//cout <<"hi";
				}
				
			}
				
		
		    else   // completely new seq_num
			{	h_index->j[hashtab->hash_val[ii]]++;	
				hv[h_index->pos[hashtab->hash_val[ii]]+h_index->j[hashtab->hash_val[ii]]-1]=ssn;
			}
	    
	   }
	   
	 }
 }
 
 	
	
		
	if(NULL != hashtab->hash_val)
		delete[] hashtab->hash_val;
	
  

	if(NULL != hashtab)
		delete hashtab; 
		
		
		
		
	
	return 1;
}


