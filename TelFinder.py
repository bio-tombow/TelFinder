
from FindKmer import *#kmer find
from statMut import *#statistic about kmers
from score import *#score
from multiprocessing import Pool

import argparse

parser = argparse.ArgumentParser(description='TelFinder is used to detect telomeric repeat sequence')
#############
parser.add_argument('-f', '--format', required=True, help='format of input genome. fa means fasta (for identifying telomere unit); kmer means kmer candidates(for mutation analysis)')
parser.add_argument('-inf','--input-file', help='input file (fa or fq genome file or kmer candidates file)')
parser.add_argument('-s', '--species', required=True, help='species for name of output file')
parser.add_argument('-r', '--size-range', default= '5,30', help='needed when -f fa or -f fq is used, size range of candidates. (Default: 5,30)')
parser.add_argument('-l', '--seq-len', default= '1000', help='needed when -f fa is used, search region (nt) from chromosome termini. (Default: = 1000)')
parser.add_argument('-m', '--motif', help='needed when -f kmer is used.')
parser.add_argument('-e', '--end', help='needed when -f kmer is used. -e five means five end, -e three means three end')
parser.add_argument('-t', '--threads',default='1',help='number of threads to launch in detection. (Default: =1)')
#############
parser.add_argument('-o', '--output-folder', default=os.getcwd(), help='Directory to write output files.  (Default: present working directory)')
args = parser.parse_args()


###########
input_file = args.input_file
input_format = args.format
species = args.species
out_dir = args.output_folder
end = args.end
###########

if input_format in ['fa','fq']:
	#################
	size_range = args.size_range
	seq_len = args.seq_len
	tnumber = args.threads
	#################
	print('input file: .fa: ')
	print('start detecting candidates of Telomeric repeat')
	#################
	size_range_info = [int(a) for a in size_range.split(',')]
	kmer_start = size_range_info[0]
	kmer_end = size_range_info[1]+1
	#################
	########input==fa#########
	if input_format == 'fa':
		chrEndSeq = GetChrEnd(int(seq_len),end,species,input_file,out_dir)
		out_file = species+'_'+end+'_info.txt'
		out1 = open(os.path.join(out_dir,out_file),'w+')
		out1.write('chr\tkmer_size\tkmer\ttotal_repeat_number\tgap_number\tgap_dis\tgap_base\n')
		for chr1 in chrEndSeq:
			#print(chr1)
			seq1 = chrEndSeq[chr1]
			component_dc, gap_dis_dc, gap_base_dc = chrKmerFind(seq1, kmer_start, kmer_end)
			#new_component_dc = filter_kmer_component(component_dc, 1)
			for kmer in component_dc:
				gap_dis_info = ''
				for each in gap_dis_dc[kmer]:
					gap_dis_info += (str(each) + ';')
				gap_base_info = ''
				for each in gap_base_dc[kmer]:
					gap_base_info += (each + ';')
				out1.write(chr1 + '\t' + str(len(kmer)) + '\t' + str(kmer) + '\t' +
					str(component_dc[kmer][0]) + '\t' + str(component_dc[kmer][1]) +
					'\t' + gap_dis_info + '\t' + gap_base_info + '\n')
		out1.close()
		chr_motif = GetMotif(os.path.join(out_dir,out_file))
		print('yesyes')
		final_motif, unimotif_ls = MotifInChr(chr_motif)
		print('yesyes!!')
		supp2 = GetSupp(unimotif_ls,final_motif)
		print('yesyes')
		print(supp2)
		out_dc1 = OutCandifa(supp2,final_motif)
		out_score = species+'_'+end+'_score.txt'
		out2 = open(os.path.join(out_dir,out_score),'w+')
		out2.write('motif\tchr_supp\tmax_score\tchrscores\n')
		for kmer in out_dc1:
			chrs,scores = out_dc1[kmer]
			#print(chrs)
			maxscore = max(scores)
			chrscores = [chrs[a]+'|'+str(scores[a]) for a in range(len(chrs))]
			chrscores = ';'.join(chrscores)
			out2.write('\t'.join([kmer,str(len(chrs))+'/'+str(len(chrEndSeq)),str(maxscore),chrscores]))
			out2.write('\n')
		out2.close()
	#######input==fq##########
	if input_format == 'fq':
		p = Pool(int(tnumber))
		###
		#os.system('hisat2-build '+)
		###
		chrEndSeq = {}
		fas = SeqIO.parse('reads.txt','fasta')##unmapped reads
		for k1 in fas:
			chrEndSeq[k1.id] = str(k1.seq)
		##
		out_file = species+'_fq_info.txt'
		out1 = open(os.path.join(out_dir,out_file),'w')
		out1.write('chr\tkmer_size\tkmer\ttotal_repeat_number\tgap_number\tgap_dis\tgap_base\n')
		out1.close()
		for chr1 in chrEndSeq:
			seq1 = chrEndSeq[chr1]
			#judge(chr1,seq1,kmer_start,kmer_end,outf)
			p.apply_async(judge,args=(chr1,seq1,kmer_start,kmer_end,outf,))
		p.close()
		p.join()
		###
		chr_motif = GetMotif(out1)
		final_motif, unimotif_ls = MotifInChr(chr_motif)
		supp2 = GetSupp(unimotif_ls,final_motif)
		out_dc1 = OutCandifq(supp2,final_motif)
		out_score = species+'_fq_score.txt'
		out2 = open(os.path.join(out_dir,out_score),'w+')
		out2.write('motif\tchr_supp\tmax_score\n')
		for kmer in out_dc1:
			out2.write(kmer+'\t'+str(out_dc1[kmer][0])+'\t'+str(out_dc1[kmer][1])+'\n')
		out2.close()
if input_format == 'kmer':
	#################
	motif = args.motif
	#################
	print('input file: kmer: '+input_file)
	print('start analysis of mutation')
	#################
	info1 = species +'_'+ end
	genomef = input_file
	chr2n,chr2len,chr2gap = IndelSeq(genomef,motif,end,out_dir)
	OutputMut(info1,motif,chr2n,chr2len,chr2gap,out_dir)


