#include <iostream>
#include <stdint.h>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include "conLSH.h"
using namespace std;


uint8_t trans[128]= {
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
		4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
		4,0,4,1,4,4,4,2,4,4,4,4,4,4,2,4,
		4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4
	};
	
	
//generate the pattern for L different hash tables
void generate_pat(uint32_t *pat, uint32_t K,uint32_t L,uint32_t lambda,uint32_t zero)
{
	uint32_t p;
	uint8_t r;
	for(uint8_t i=0;i<L;i++)
	{
		pat[i]=0;
		p=(1<<lambda)-1;
		uint8_t kcount=1;
		while(kcount<=K)
		{
			pat[i]=(pat[i]<<lambda)|p;
			if(i==0)
				 r=(rand()%(5));
			else
				 r=(rand()%(zero+1));
			pat[i]=pat[i]<< r;
			kcount++;
		}
		//cout<<pat[i]<<endl;;
	}
}


// Generates the conLSH values for the entire Genome and returns the no. hash values generated
int generate_conLSH( uint32_t *h_pos, uint8_t *h_j, uint32_t *hval, uint32_t len_genome, uint32_t len_sed, char *genome, uint32_t L,uint32_t *pat)
{
	uint32_t  acn=0;
	uint32_t num_pre;
	uint32_t  cc=1;  // cc stores count of valid sequences excluding the N's and overrepresented ones
	uint32_t j,i;
	
	for( i = 1; i < len_genome-len_sed; i+=1) {
		
		
		
		//skip the N's
		
		for (j=0;j+i<len_genome-len_sed && genome[i+j]=='N';++j);
		if (i+j > len_genome - len_sed) break;
		i+=j;
		
		//if(ww>50) //window size needs to be a parameter
		//{
			//ww=1;
			//seq++;
		//}
	
	    // Use only one pattern for genome
		
			transfer(genome + i,&num_pre,len_sed,pat[0]);
			
			
			if(h_j[num_pre]<=199)
			{
				//cout<<num_pre<<"\t"<<h_j[num_pre]<<endl;
				hval[i]=num_pre;
				//sn[count]=i;
			
			  
			   if(h_j[num_pre]==0)
			   {
				   h_pos[num_pre]=0;	
					cc++;
		       }
		        h_j[num_pre]++;
		       acn++;
		      
			}
			else  //over-represented hash_value
			{
				//cout<<num_pre<<"\t"<<h_j[num_pre]<<endl;
				if(h_j[num_pre]==200)
				{
					acn=acn-h_j[num_pre];
					h_j[num_pre]=201;        
				}
				
			}
			//cout<<"h_index="<<h_index[num_pre]<<endl;
			//if(num_pre>mx)
				//mx=num_pre;
			//cout <<seq<<"\t"<<num_pre<<"\n";
			
			
		}
		
		
	
	
	//*c=acn-1;
	//cout<<"distinct hash_value="<<cc-1<<"\t total valid sequences="<<acn-1<<"\n";
	return(acn);
}


// Generates the conLSH values for the long and noisy reads
int generate_conLSH_read(uint32_t *hv, uint32_t len_read, char *read, uint32_t L,uint32_t *pat, int tr)
{
	uint32_t count=0;
	uint32_t num_pre, amount;
	
	if(tr==4)
	amount=((len_read - 32)>>2)*L ;
	else if(tr==8)
	amount=((len_read - 32)>>3)*L ;
	else if (tr==16)
	 amount=((len_read - 32)>>4)*L ;
	 else if (tr==32)
	 amount=((len_read - 32)>>5)*L ;
	 else if (tr==64)
	 amount=((len_read - 32)>>6)*L ;
	 else if (tr==128)
	 amount=((len_read - 32)>>7)*L ;
	 
	for (uint32_t pp=0;pp<L;pp++) // generate hash values for L different patterns
   {
		
	for( uint32_t i = 1; i < len_read - 32; i+=tr) {
		
		
		uint32_t j;
		//skip the N's
		
		for (j=0;j+i<len_read-32 && read[i+j]=='N';++j);
		if (i+j > len_read - 32) break;
		i+=j;
		
		
		
			
			transfer(read + i,&num_pre,32,pat[pp]);
			if(count<amount)
			{
			  hv[count]=num_pre;
			count++;
			}
		}
		
		
	}
	return(count);
}

int transfer(char *genome,uint32_t *l2r,uint32_t len_sed, uint32_t pt)
{

	*l2r = 0;
	uint32_t m;
	m=1<<(len_sed-1);
	uint32_t temp;
	for(uint32_t i = 0;i<len_sed;i++) {
		
		if(pt&m)
		{
		temp = trans[genome[i]];
		*l2r = *l2r << 2;
		*l2r = *l2r | temp;
		}
		m=m>>1;
	}
	return 0;
}



