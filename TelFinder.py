
from FindKmer import *#kmer find
from statMut import *#statistic about kmers
from score import *#score
from multiprocessing import Pool

import argparse

parser = argparse.ArgumentParser(description='TelFinder is used to detect telomeric repeat sequence')
#############
parser.add_argument('-f', '--format', required=True, help='format of input genome. fa means fasta, fq1 means single end fastq, fq2 means paired end fastq (for identifying telomere unit); kmer means kmer candidates(for mutation analysis)')
parser.add_argument('-inf','--input-file',help='input file (fa or fq(single end) genome file or kmer candidates file)')
parser.add_argument('-s', '--species', required=True, help='species for name of output file')
parser.add_argument('-r', '--size-range', default= '5,30', help='size range of candidates. (Default: 5,30)')
parser.add_argument('-l', '--seq-len', default= '1000', help='needed when -f fa is used, search region (nt) from chromosome termini. (Default: = 1000)')
parser.add_argument('-m', '--motif', help='needed when -f kmer is used.')
parser.add_argument('-e', '--end', default = 'left',help='-e left or -e right can be set. (Default: -e left)')
parser.add_argument('-o', '--output-folder', default=os.getcwd(), help='Directory to write output files.  (Default: current working directory)')
###
parser.add_argument('-1', '--fq-p1', help='needed when -f fq2 is used')
parser.add_argument('-2', '--fq-p2', help='needed when -f fq2 is used')
parser.add_argument('-ref', '--reference', help='needed when -f fq1/fq2 is used')
###
parser.add_argument('-t', '--threads',default='1',help='number of threads to launch in detection. (Default: =1)')
#############
args = parser.parse_args()


###########
input_format = args.format
species = args.species
out_dir = args.output_folder
tnumber = args.threads
###########

if input_format in ['fa','fq1','fq2']:
	########input==fa#########
	if input_format == 'fa':
		#################
		input_file = args.input_file
		size_range = args.size_range
		seq_len = args.seq_len
		end = args.end
		#################
		print('input file: .fa: ')
		print('start detecting candidates of Telomeric repeat')
		#################
		size_range_info = [int(a) for a in size_range.split(',')]
		kmer_start = size_range_info[0]
		kmer_end = size_range_info[1]+1
		#################
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
		final_motif, unimotif_ls = MotifInChr(chr_motif)
		supp2 = GetSupp(unimotif_ls,final_motif)
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
	if 'fq' in input_format:
		size_range = args.size_range
		size_range_info = [int(a) for a in size_range.split(',')]
		kmer_start = size_range_info[0]
		kmer_end = size_range_info[1]+1
		tnumber = args.threads
		p = Pool(int(tnumber))
		#########
		ref= args.reference
		print('build reference...')
		os.system('hisat2-build -q '+ref+' '+os.path.join(out_dir,species))
		print('alignment...')
		if input_format == 'fq1':
			input_file = args.input_file
			mode = 'single'
			os.system('hisat2 -p '+tnumber +' -x '+os.path.join(out_dir,species)+' -U '+input_file+' -S '+os.path.join(out_dir,species+'.sam'))
		elif input_format == 'fq2':
			mode = 'paired'
			p2f1 = args.fq_p1
			p2f2 = args.fq_p2
			os.system('hisat2 -p '+tnumber +' -x '+os.path.join(out_dir,species)+' -1 '+p2f1+' -2 '+p2f2+' -S '+os.path.join(out_dir,species+'.sam'))
		#########
		print('extract unmapreads...')
		exunmapped(os.path.join(out_dir,species+'.sam'),mode,os.path.join(out_dir,species+'_unmapreads.fasta'))
		#########
		chrEndSeq = {}
		fas = SeqIO.parse(os.path.join(out_dir,species+'_unmapreads.fasta'),'fasta')##unmapped reads
		for k1 in fas:
			chrEndSeq[k1.id] = str(k1.seq)
		##
		print('Motif search...')
		out_file = species+'_fq_info.txt'
		out1 = open(os.path.join(out_dir,out_file),'w')
		out1.write('chr\tkmer_size\tkmer\ttotal_repeat_number\tgap_number\tgap_dis\tgap_base\n')
		out1.close()
		for chr1 in chrEndSeq:
			seq1 = chrEndSeq[chr1]
			#judge(chr1,seq1,kmer_start,kmer_end,os.path.join(out_dir,out_file))
			p.apply_async(judge,args=(chr1,seq1,kmer_start,kmer_end,os.path.join(out_dir,out_file),))
		p.close()
		p.join()
		###
		chr_motif = GetMotif(os.path.join(out_dir,out_file))
		final_motif, unimotif_ls = MotifInChr(chr_motif)
		supp2 = GetSupp(unimotif_ls,final_motif)
		out_dc1 = OutCandifq(supp2,final_motif)
		out_score = species+'_fq_score.txt'
		out2 = open(os.path.join(out_dir,out_score),'w+')
		out2.write('motif\tchr_supp\tmax_score\n')
		for kmer in out_dc1:
			out2.write(kmer+'\t'+str(out_dc1[kmer][0])+'\t'+str(out_dc1[kmer][1])+'\n')
		out2.close()
if 'kmer' in input_format:
	#################
	end = args.end
	motif = args.motif
	if input_format == 'kmer_fa':
		input_file = args.input_file
		#################
		print('input file: kmer: '+input_file)
		print('start analysis of mutation')
		#################
		info1 = species +'_'+ end
		genomef = input_file
		chr2n,chr2len,chr2gap = IndelSeq(genomef,motif,end,out_dir)
		OutputMut(info1,motif,chr2n,chr2len,chr2gap,out_dir)
	if 'kmer_fq' in input_format:
		###
		ref= args.reference
		print('build reference...')
		os.system('hisat2-build -q '+ref+' '+os.path.join(out_dir,species))
		print('alignment...')
		if input_format[-1] == '1':
			input_file = args.input_file
			print('input file: kmer: '+input_file)
			mode = 'single'
			os.system('hisat2 -p '+tnumber +' -x '+os.path.join(out_dir,species)+' -U '+input_file+' -S '+os.path.join(out_dir,species+'.sam'))
		elif input_format[-1] == '2':
			mode = 'paired'
			p2f1 = args.fq_p1
			p2f2 = args.fq_p2
			print('input file: kmer: '+','.join([p2f1,p2f2]))
			os.system('hisat2 -p '+tnumber +' -x '+os.path.join(out_dir,species)+' -1 '+p2f1+' -2 '+p2f2+' -S '+os.path.join(out_dir,species+'.sam'))
		#########
		print('extract unmapreads...')
		unmapf = os.path.join(out_dir,species+'_unmapreads.fasta')
		exunmapped(os.path.join(out_dir,species+'.sam'),mode,unmapf)
		###
		print('start analysis of mutation')
		telreadsf = os.path.join(out_dir,'telreads')
		h1 = open(telreadsf,'w')
		with open(unmapf) as f:
			for lines in f.readlines():
				if lines[0] == '>':
					id1 = lines[1:-1]
				else:
					seq1 = lines
					if motif*3 in seq1:
						h1.write('>'+id1+'\n'+seq1)
		h1.close()
		os.system('rm '+unmapf)
		###
		info1 = species+'_'+end
		chr2n,chr2len,chr2gap = IndelSeq(telreadsf,motif,end,out_dir)
		OutputMutfq(info1,motif,chr2n,chr2len,chr2gap,out_dir)
		os.system('rm '+telreadsf)


