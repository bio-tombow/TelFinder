------------------------------------------------------
TelFinder: a program for de novo detection of telomeric motif sequence
python3 is needed
If you have any problem in using the program, please no hesitate to 
contact Dr.Sun: tombow.bioinfo@gmail.com
------------------------------------------------------

Usage:
cd TelFinder-main/

#check the parameters and its meanings
python3 TelFinder -h

#format=fasta,
python3 TelFinder -f fa -inf example/fasta/NCR_genomic.fna -s NCR -o example/fasta -e five
python3 TelFinder -f fa -inf example/fasta/NCR_genomic.fna -s NCR -o example/fasta -e three

#format=fastq,
python3 TelFinder -f fq -inf example/fastq/test.genomic.fastq -s test -o example/fastq

#varaition analysis, format=kmer
python3 TelFinder -f kmer -inf example/fasta/ara_five_info -m TTAGGG -o example/fasta -e five

