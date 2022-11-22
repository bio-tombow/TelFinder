
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
            times = int(line[3])
            if line[4] !='0':
                for gap in line[-1][:-1].split(';'):
                    gap_n=len(gap)
                    #print(gap)
                    if gap_n >= len(motif)*3 and motif not in gap:
                        times -=1
            t1 = [motif, float(line[1])*1, float(times)*1]#motif,length,times
            if len(line) == 7:
                t1.append(line[6])
            else:
                t1.append('')
            if line[0] not in chr_motif:
                chr_motif[line[0]] = [t1]
            else:
                chr_motif[line[0]].append(t1)
    return chr_motif

def minseq(my_seq):
    minseq = my_seq
    for i in range(1,len(my_seq)):
        if len(my_seq) % i ==0 and i !=1:
            #print(i)
            tmp = []
            for s1 in range(0,len(my_seq),i):
                seq1 = my_seq[s1:s1+i]
                #print(seq1)
                if seq1 not in tmp:
                    tmp.append(seq1)
            if len(tmp) == 1:
                minseq = tmp[0]
                #print('minseq:'+minseq)
                break
    return minseq

##########
def ifsub(motif1,motif2):#一个motif是否是另一个motif的整倍数
    out = 'None'
    if len(motif1) > len(motif2):
        a1 = len(motif1) / len(motif2)
        motif_2 = motif2*int(a1)
        if compare_kmer(motif1,motif_2) =='same':
            out = motif2
    else:
        a1 = len(motif2) / len(motif1)
        motif_1 = motif1*int(a1)
        if compare_kmer(motif_1,motif2) =='same':
            out = motif1
    return out

##########

def MotifInChr(chr_motif):
    final_motif = {}
    unimotif_ls = []
    for chr1 in chr_motif:
        motif_ls = chr_motif[chr1]
        print('build motif2score')
        motif2score={}
        tmp_ls1 = []
        supp1 = {}
        for each in motif_ls:
            supp1[each[0]] = 0
        if len(chr_motif[chr1]) ==1:
            motif_info = chr_motif[chr1][0]
            motif2score[motif_info[0]] = motif_info[1] * motif_info[2]
            supp1[motif_info[0]] = 1
        else:
            for motif_info1 in chr_motif[chr1]:
                tmp_ls1.append(motif_info1)
                for motif_info2 in motif_ls:
                    print(motif_info2,'motifinfo2')
                    if motif_info2 not in tmp_ls1:
                        motif1 = motif_info1[0]
                        score1 = motif_info1[1] * motif_info1[2]  #
                        motif2 = motif_info2[0]
                        score2 = motif_info2[1] * motif_info2[2]  #
                        print(motif1,motif_info1,'tmp')
                        #print(motif_info2)
                        #print(motif1 + '-' + motif2)
                        #print(score1)
                        #print(score2)
                        score = score1 - score2
                        if motif_info1[-1] != "":  # 为消除子串影响
                            tmp1 = 0
                            gaps = motif_info1[-1].split(';')[:-1]
                            for gap1 in gaps:
                                comb = motif1+gap1
                                if compare_kmer(comb,motif2)!= 'diff':#相同或为子串
                                    tmp1 +=1
                                    #print(motif1+'-'+motif2)
                                    #print(gap1)
                            score = score - tmp1 * motif_info1[1]#len * 系数
                        #print(motif1 + '-' + motif2 + '|' + str(score))
                        #print(score)
                        ######################
                        if motif1 not in motif2score:
                            motif2score[motif1] = score1
                        else:
                            if score1 >= motif2score[motif1]:
                                motif2score[motif1] = score1
                        if motif2 not in motif2score:
                            motif2score[motif2] = score2
                        else:
                            if score2 >= motif2score[motif2]:
                                motif2score[motif2] = score2
                        #######################
                        out = 'None'
                        if score1 == score2:
                            out = ifsub(motif1,motif2)
                        if out != 'None':
                            supp1[out] +=1
                        else:
                            if score > 0:
                                supp1[motif1] += 1
                            else:
                                supp1[motif2] += 1
        print(chr1)
        print('motif2score end')
        print(motif2score)
        print(supp1)
        ##################
        del1 = []
        mins = []
        for k1 in supp1:
            if minseq(k1) ==k1:
                mins.append(k1)
            else:
                for min1 in mins:
                    if compare_kmer(min1,minseq(k1)) == 'same':
                        del1.append(k1)
        for k2 in del1:
            try:
                del supp1[k2]
            except:
                pass
        print(supp1)
        print('\n')
        #################
        tmp_ls2 = [a for a in supp1.values()]
        for k3 in supp1:
            if max(tmp_ls2) ==0:
                pass
            else:
                if supp1[k3] == max(tmp_ls2):
                    if k3 not in unimotif_ls:
                        flag = 0
                        for k4 in unimotif_ls:
                            if compare_kmer(k3, k4) == 'same':
                                flag = 1
                        if flag == 0:
                            unimotif_ls.append(k3)
                    #print(k3)
                    #print(motif2score[k3])
                    if chr1 not in final_motif:
                        final_motif[chr1] = [[k3,motif2score[k3]]]
                    else:
                        final_motif[chr1].append([k3,motif2score[k3]])
    return final_motif, unimotif_ls

def GetSupp(unimotif_ls,final_motif):
    #print(unimotif_ls)
    #print(final_motif)
    supp2 = {}
    for each in unimotif_ls:
        print(each)
        supp2[each] = []
        for chr1 in final_motif:
            flag = 0
            for motif,score in final_motif[chr1]:
                if compare_kmer(each, motif) == 'same':
                    flag = 1
            if flag == 1:
                #print(each)
                print(chr1)
                #supp2[each] += 1
                supp2[each].append(chr1)
    return supp2


def OutCandifa(supp2,final_motif):
    ######
    kmer2cs = {}#kmer 2 chrn,maxscore
    kmer2chrs = {}
    kmer2scores = {}
    for each in supp2:
        kmer2cs[each] = [0,0]
        kmer2cs[each][0] = len(supp2[each])
        kmer2chrs[each] = supp2[each]
    for each in kmer2cs:
        tmp_scores = []
        for chr1 in final_motif:
            for motif,score in final_motif[chr1]:
                if compare_kmer(each,motif) == 'same':
                    tmp_scores.append(round(score,3))
        kmer2scores[each] = tmp_scores
        kmer2cs[each][1] = max(tmp_scores)
    #######
    out_dc1 = {}
    kmer2cs1 = sorted(kmer2cs.items(), key=lambda x: x[1], reverse=True)
    for each in kmer2cs1:
        kmer = each[0]
        chrs = kmer2chrs[kmer]
        scores = kmer2scores[kmer]
        out_dc1[kmer] = [chrs,scores]
    #######
    return out_dc1


def OutCandifq(supp2,final_motif):
    #print(supp2)
    values = [a for a in supp2.values()]
    values.sort(reverse=True)
    ##############
    kmer2cs = {}#kmer 2 chrn,maxscore
    for each in supp2:
        kmer2cs[each] = [0,0]
        kmer2cs[each][0] = len(supp2[each])
    for each in kmer2cs:
        tmp_scores = []
        for chr1 in final_motif:
            for motif,score in final_motif[chr1]:
                if compare_kmer(each,motif) == 'same':
                    tmp_scores.append(round(score,3))
        kmer2cs[each][1] = max(tmp_scores)
    #######
    out_dc1 = {}
    kmer2cs1 = sorted(kmer2cs.items(), key=lambda x: x[1], reverse=True)
    for each in kmer2cs1:
        kmer = each[0]
        chrn,maxscore = kmer2cs[kmer]
        out_dc1[kmer] = [chrn,maxscore]
    return out_dc1
