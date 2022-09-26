import collections
from Bio import SeqIO
import os
#########################


def GetChrEnd(len0,end,species,fasta_file,dir1):
    chrEndSeq = {}
    which = ['five', 'three']
    index = which.index(end)
    info1 = []
    fas = SeqIO.parse(fasta_file, 'fasta')
    for k1 in fas:
        id1 = k1.id
        seq = str(k1.seq).upper()
        seq = seq.strip('N')
        if info1 == []:
            info1.append(id1[:2])
        if id1[:2] == info1[0]:
            if index ==0:
                chrEndSeq[k1.id] = seq[:len0].upper()
            if index ==1:
                chrEndSeq[k1.id] = seq[-(len0):].upper()
    h1 = open(os.path.join(dir1,'_'.join([species,end])),'w+')
    for k1 in chrEndSeq:
        h1.write('>'+k1+'\n'+chrEndSeq[k1]+'\n')
    h1.close()
    return chrEndSeq

def split_dna(dna, kmer_size):
    kmers = []
    for start1 in range(0,len(dna)-(kmer_size-1)):
        kmer1 = dna[start1:start1+kmer_size]
        kmers.append(kmer1)
    return kmers

#########################
###if kmer1 == TTAGGG,kmer2 == GGGTTA, then out='same'.compare 2 kmers
def compare_kmer(kmer1,kmer2):
	out = ''
	if len(kmer1) == len(kmer2):
		new_kmer = kmer1+kmer1
		if kmer2 == kmer1 or kmer2 in new_kmer:
			out = 'same'
		else:
			out = 'diff'
	else:
		if len(kmer1) > len(kmer2):
			long_kmer = kmer1
			short_kmer = kmer2
		else:
			long_kmer = kmer2
			short_kmer = kmer1
		new_long_kmer = long_kmer+long_kmer
		if short_kmer in long_kmer or short_kmer in new_long_kmer:
			out = 'similar'
		else:
			out = 'diff'
	return out

#########################
## reverse_transcript, emm, no reverse
def RT(seq1):
	seq2 = ''
	dc = {'A':'T','T':'A','C':'G','G':'C'}
	for k in seq1:
		seq2 += dc[k]
	return seq2

#########################
## remove kmer only with base of one type: TTTT/AAAA/CCCC/GGGGG
def rm1base(kmer_list):
    site_list1 = []
    kmer_list1 = []
    kmer_dict1 = {}
    for site in range(len(kmer_list)):
        each = kmer_list[site]
        each1 = list(each)
        counter = collections.Counter(each1)
        if len(counter) >2:
            kmer_list1.append(each)
            site_list1.append(site)
            kmer_dict1[site] = each

    return site_list1, kmer_list1, kmer_dict1

#########################   
## count the number of non-redundant repeat, 
#like TTAGGG-T-TTAGGG-TTAGGG-G-TTAGGG, result will be TTAGGG:4
# new_kmers[TTAGGG] = 4; new_sites[TTAGGG] = [0,7,13,]
def count_kmer(kmer_list, kmer_dict):
    count1 = collections.Counter(kmer_list)
    kmer_repeat = {}
    kmer_site = {}
    kmer_site1 = {}
    for i1 in count1:
        if count1[i1] != 1:
            match = []
            for i2 in kmer_repeat:
                out = compare_kmer(i1, i2)
                if out == 'same':
                    if count1[i1] > kmer_repeat[i2]:
                        del kmer_repeat[i2]
                        kmer_repeat[i1] = count1[i1]
                        break
                match.append(out)
            if 'same' not in match:
                # print('add ' + k1)
                kmer_repeat[i1] = count1[i1]
    for i1 in kmer_repeat:
        kmer_site[i1] = []
        kmer_site1[i1] = []
        for i2 in kmer_dict:
            if kmer_dict[i2] == i1:
                kmer_site[i1].append(i2)
                kmer_site1[i1].append(i2)
    for i1 in kmer_site:
        ls1 = kmer_site[i1]
        ls2 = kmer_site1[i1]
        for i2 in range(len(ls1)):
            for i3 in range(len(ls2)):
                if i3 >0 and i3 < len(ls2):
                    if ls2[i3-1] + len(i1) > ls2[i3]:
                        ls2.remove(ls2[i3])
                        break
        kmer_site[i1] = ls2
        kmer_repeat[i1] = len(ls2)
    return kmer_repeat, kmer_site

#########################
def kmer_info(kmer_size,kmer_site,kmer_repeat,component_dc,gap_dis_dc,gap_base_dc,seq1):
    for k1 in kmer_repeat:
        component_dc[k1] = [kmer_repeat[k1], 0]
        #repeat time, gap number
        gap_dis_dc[k1] = []
        #gap dis list
        gap_base_dc[k1] = []
        #gap base list
        ls = [kmer_site[k1][0]]
        for k2 in kmer_site[k1][1:]:
            if k2 > ls[-1] + kmer_size:
                gap_base = seq1[ls[-1]+kmer_size:k2]
                gap_dis = len(gap_base)
                component_dc[k1][1] += 1
                gap_dis_dc[k1].append(gap_dis)
                gap_base_dc[k1].append(gap_base)
            if k2 == ls[-1] + kmer_size:
                pass
            ls.append(k2)
        if component_dc[k1][0] == component_dc[k1][1]+ 1:
            del component_dc[k1]
            del gap_dis_dc[k1]
            del gap_base_dc[k1]
    return component_dc, gap_dis_dc, gap_base_dc

#########################
def filter_kmer_component(component_dc, repeat_ratio):
    repeat_time1 = []
    for k1 in component_dc:
        repeat_time1.append(component_dc[k1][0])
    repeat_time2 = list(set(repeat_time1))
    repeat_time2.sort(reverse=True)
    threshold = repeat_time2[:int(len(repeat_time2) * repeat_ratio)]
    new_component_dc = {}
    for k1 in component_dc:
        if threshold == []:
            pass
        else:
            if component_dc[k1][0] > min(threshold):
                new_component_dc[k1] = component_dc[k1]
    return new_component_dc

#########################
def chrKmerFind(seq1, kmer_start, kmer_end):
    component_dc = {}
    gap_dis_dc = {}
    gap_base_dc = {}
    for kmer_size in range(kmer_start, kmer_end):
        #print(kmer_size)
        kmer_list = split_dna(seq1, kmer_size)
        site_list1, kmer_list1, kmer_dict1 = rm1base(kmer_list)
        kmer_repeat, kmer_site = count_kmer(kmer_list1, kmer_dict1)
        component_dc, gap_dis_dc, gap_base_dc = kmer_info(kmer_size, kmer_site, kmer_repeat,
                                                          component_dc, gap_dis_dc, gap_base_dc, seq1)
    return component_dc, gap_dis_dc, gap_base_dc

#########################

