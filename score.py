
import os

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

def GetMotif(info_file):
    chr_motif = {}
    with open(info_file) as f1:
        for lines in f1.readlines()[1:]:
            line = lines.strip().split('\t')
            motif = line[2]
            if len(line) == 7:
                t1 = [motif, line[1], line[3], line[6]]
            else:
                t1 = [motif, line[1], line[3], '']
            if line[0] not in chr_motif:
                chr_motif[line[0]] = [t1]
            else:
                chr_motif[line[0]].append(t1)
    return chr_motif

def MotifInChr(chr_motif):
    final_motif = {}
    unimotif_ls = []
    for chr in chr_motif:
        motif_ls = chr_motif[chr]
        tmp_ls1 = []
        supp1 = {}
        for each in motif_ls:
            supp1[each[0]] = 0
        for motif_info1 in chr_motif[chr]:
            tmp_ls1.append(motif_info1)
            for motif_info2 in motif_ls:
                if motif_info2 not in tmp_ls1:
                    motif1 = motif_info1[0]
                    score1 = int(motif_info1[1]) * int(motif_info1[2])  #
                    motif2 = motif_info2[0]
                    score2 = int(motif_info2[1]) * int(motif_info2[2])  #
                    #print(motif1 + '-' + motif2)
                    #print(score1)
                    #print(score2)
                    score = score1 - score2
                    if motif_info1[-1] != "":  # 为消除子串影响
                        tmp1 = 0
                        gaps = motif_info1[-1].split(';')[:-1]
                        for gap1 in gaps:
                            if motif1 + gap1 == motif2:
                                #print(gap1)
                                tmp1 += 1
                            elif len(motif1) + len(gap1) >= len(motif2):
                                tmp2 = len(motif2) - len(motif1)
                                if motif1 + gap1[:tmp2] == motif2:
                                    #print(gap1)
                                    tmp1 += 1
                        score = score - tmp1 * len(motif1)
                    #print(motif1 + '-' + motif2 + '|' + str(score))
                    #print(score)
                    if score > 0:
                        supp1[motif1] += 1
                    else:
                        supp1[motif2] += 1
        del1 = []
        for k1 in supp1:
            for i in range(2, 6):
                if k1 * i in supp1:
                    del1.append(k1 * i)
        for k2 in del1:
            try:
                del supp1[k2]
            except:
                print(supp1)
                print(k2)
        tmp_ls2 = [a for a in supp1.values()]
        for k3 in supp1:
            if supp1[k3] == max(tmp_ls2):
                if k3 not in unimotif_ls:
                    flag = 0
                    for k4 in unimotif_ls:
                        if compare_kmer(k3, k4) == 'same':
                            flag = 1
                    if flag == 0:
                        unimotif_ls.append(k3)
                if chr not in final_motif:
                    final_motif[chr] = [k3]
                else:
                    final_motif[chr].append(k3)
    return final_motif, unimotif_ls

def GetSupp(unimotif_ls,final_motif):
    supp2 = {}
    for each in unimotif_ls:
        supp2[each] = 0
        for chr in final_motif:
            flag = 0
            for motif in final_motif[chr]:
                if compare_kmer(each, motif) == 'same':
                    flag = 1
            if flag == 1:
                supp2[each] += 1
    return supp2


def OutCandi(supp2,chr_motif):
    values = [a for a in supp2.values()]
    values.sort(reverse=True)
    out_dc1 = {}
    for each in supp2:
        if supp2[each] >= values[1]:
            out_dc1[each] = [supp2[each]]
    for each in out_dc1:
        tmp_scores = []
        for chr in chr_motif:
            for motif_info1 in chr_motif[chr]:
                if compare_kmer(each, motif_info1[0]) == 'same':
                    s1 = int(motif_info1[1]) * int(motif_info1[2])
                    tmp_scores.append(s1)
        out_dc1[each].append(max(tmp_scores))
    return out_dc1

