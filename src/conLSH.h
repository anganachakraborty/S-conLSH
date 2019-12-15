


//generate the pattern for L different hash tables
void generate_pat(uint32_t *pat, uint32_t K,uint32_t L,uint32_t lambda,uint32_t zero);

// Generates the conLSH values for the entire Genome and returns the no. hash values generated
int generate_conLSH( uint32_t *h_pos, uint8_t *h_j, uint32_t *hval, uint32_t len_genome, uint32_t len_sed, char *genome, uint32_t L,uint32_t *pat);

// Generates the conLSH values for the short reads
int generate_conLSH_read(uint32_t *hv, uint32_t len_read, char *read, uint32_t L,uint32_t *pat, int tr);

int transfer(char *genome,uint32_t *l2r,uint32_t len_sed, uint32_t pt);



