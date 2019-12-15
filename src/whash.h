#ifndef WHASH_H_
#define WHASH_H_
	#define NUM_FILE 256
	#define N_LIMIT 9
	#define LEN_BASES 1024

	#define LEN 100

	#include <stdint.h>
        #include "conLSH.h" 
	typedef struct{
		//uint32_t *sqn; //seq_number
		uint32_t *hash_val; // conLSH value
	}Hashtab;
        typedef struct{
		uint32_t *pos; //start position of seq in that hash_value (initially the count of distict hash values)
		uint8_t *j; // last inserted position in hv (initially the frequency of that hash value)
	}H_index;

	
	class Hash{
	public:
                    uint32_t *pat;	
		    Hashtab * hashtab;
		     H_index * h_index;
		     uint32_t *hv;
			uint32_t count;
		int write_hashfile(char *path,char *genome,uint32_t len_genome,uint32_t len_sed,uint32_t K,uint32_t L,uint32_t lambda,uint32_t zero);
                
	
	};
	extern uint8_t trans[];
	int transfer(char *genome,uint32_t *l2r,uint32_t len_sed, uint32_t pt);
#endif
