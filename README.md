# S-conLSH
Spaced-context based Locality Sensitive Hashing

#ntroduction
```

Spaced context based Locality Sensitive Hashing (S-conLSH) is a new mapper that facilitates gapped mapping of noisy long reads to the corresponding target locations of a reference genome, With multiple spaced patterns. We have examined the performance of the proposed method on 5 different real and simulated datasets.
S-conLSH is at least 2 times faster than the state-of-the-art alignment-based methods. It achieves a sensitivity of 99%, without using any traditional base-to-base alignment, on human simulated sequence data. By default, S-conLSH provides an alignment-free mapping in PAF format. If a base level alignment is required, S-conLSH provides an option (--align 1) to generate alignment in SAM format using ksw library (https://github.com/attractivechaos/klib).


S-conLSH is open source and free for non-commercial use.

S-conLSH is designed by Angana Chakraborty in collaboration with Sanghamitra Bandyopadhyay, Indian Statistical Institute, Kolkata and Prof. Burkhard Morgenstern, University of Göttingen, Germany. 

---
```


###Installation
```

Current version of S-conLSH needs to be run on Linux operating system.

The source code is written in C++. 

The makefile is attached. Use the make command for generating the executable file.
The binary 'S-conLSH' performs indexing of the reference genome and then aligns the long and noisy PacBio reads to it.

---
```
###Synopsis
```

S-conLSH <PathOfSourceFiles> <ReferenceGenome> <ReadFile>  [-K concatenationFactor] [-L NumberOfHashTables] [--lambda contextFactor] [--zero spacesInPatterns] [-w windowsHits] [-m candidates] [-x match] [-y mismatch] [-q gapOpen] [-r gapExtension] [-a alignInSAM] > <OutputFile>


---
```
###Parameters (could be updated in the future for adding new functions)
```
------------------------------------------------------------------------------------------------------
-K, --K                <int>           Concatenation factor of locality sensitive hashing 
-L, --L                <int>           Number of hash tables in conLSH framework
--lambda               <int>           The context factor
--zero                 <int>           The number of don't cares or zeros in the S-conLSH pattern 
-w, --window-hits      <int>           The max allowed number of windows hitting by a k-mer [Default=1000] 
-m, --candidates       <int>           The number of candidates for extension [Default=400]
-x, --match            <int>           Score of match for the alignments in extension phase [Default=2]
-y, --mismatch         <int>           Mismatch penalty for the alignments in extension phase [Default=5]
-q, --gap-open         <int>           Gap open penalty for the alignments in extension phase [Default=2]
-r, --gap-extension    <int>           Gap extension penalty for the alignments in extension phase [Default=1]
-a, --align	       <int>           Value=1, outputs alignment in SAM format [Default=0, Alignment Free PAF format output]
-h, --help                             Help
-------------------------------------------------------------------------------------------------------



---

###Quick start
```
The package includes sample reference genome and SMRT reads to demonstrate the quich start guide. 

``` For alignment free mapping in PAF format

./S-conLSH ../src/ ../sample_data/ecoli_AE005174v2.fas ../sample_data/SRR801638.fasta > sample.paf

``` For SMRT aignment in SAM format

``` ./S-conLSH ../src/ ../sample_data/ecoli_AE005174v2.fas ../sample_data/SRR801638.fasta --align 1 > sample.sam

---


###Reference

rHAT: [1]Liu, B., Guan, D., Teng, M. & Wang, Y. rHAT: fast alignment of noisy long reads with regional hashing. 
Bioinformatics 32, 1625–1631 (2015).
conLSH: [2]conLSH:Context based Locality Sensitive Hashing for Mapping of noisy SMRT Reads. Computational Biology and Chemistry, Elsevier [Accepted]

---

###Contact

For advising, bug reporting and requiring help, please contact angana_r@isical.ac.in
