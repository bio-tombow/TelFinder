import os
import math
from Bio import SeqIO

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

def gapfind(seq1,motif,end):
    motif_site = []
    gaps2seq = {}
    if end == 'left':
        i=0
        while i < len(seq1):
            if motif==seq1[i:i+len(motif)]:
                motif_site.append([i,i+len(motif)-1])
                if len(motif_site) == 1:
                    pass
                    #h1.write(','.join(['motif',str(i),str(i+len(motif)-1),str(len(motif))])+'\t'+seq1[i:i+len(motif)]+'\n')
                if len(motif_site) >=2:
                    if motif_site[-1][0] - motif_site[-2][1] <1000:
                        #h1.write(','.join(['motif',str(i),str(i+len(motif)-1),str(len(motif))])+'\t'+seq1[i:i+len(motif)]+'\n')
                        if motif_site[-1][0] - motif_site[-2][1] >1:
                            gs1 = motif_site[-2][1]+1
                            ge1 = motif_site[-1][0]-1
                            #print(gs1,ge1)
                            gaps2seq[','.join([str(gs1),str(ge1)])] = seq1[gs1:ge1+1]
                            #h1.write(','.join(['gap',str(gs1),str(ge1),str(ge1-gs1+1)])+'\t'+seq1[gs1:ge1+1]+'\n')
                    else:
                        break
                i = i + len(motif)
            else:
                i+=1
    if end == 'right':
        i=-1
        #h1 = open('three.txt','w')
        while i > 0-len(seq1):
            if motif==seq1[i-len(motif):i]:
                motif_site.append([i-len(motif)+1,i])
                if len(motif_site) ==1:
                    pass
                    #h1.write(','.join(['motif',str(0-i),str(0-(i-len(motif)+1)),str(len(motif))])+'\t'+seq1[i-len(motif):i]+'\n')
                    #h1.write(','.join(['motif',str(i-len(motif)+1),str(i),str(len(motif))])+'\t'+seq1[i-len(motif):i]+'\n')
                if len(motif_site) >=2:
                    if abs(motif_site[-1][1]) - abs(motif_site[-2][0]) <1000:
                        #h1.write(','.join(['motif',str(0-i),str(0-(i-len(motif)+1)),str(len(motif))])+'\t'+seq1[i-len(motif):i]+'\n')
                        #h1.write(','.join(['motif',str(i-len(motif)+1),str(i),str(len(motif))])+'\t'+seq1[i-len(motif):i]+'\n')
                        if abs(motif_site[-1][1]) - abs(motif_site[-2][0]) >1:
                            gs1 = motif_site[-1][1]+1
                            ge1 = motif_site[-2][0]-1
                            #print(gs1,ge1)
                            gaps2seq[','.join([str(0-ge1),str(0-gs1)])] = seq1[gs1:ge1+1]
                            #h1.write(','.join(['gap',str(0-ge1),str(0-gs1),str(ge1-gs1+1)])+'\t'+seq1[gs1:ge1+1]+'\n')
                            #h1.write(','.join(['gap',str(gs1),str(ge1),str(ge1-gs1+1)])+'\t'+seq1[gs1:ge1+1]+'\n')
                    else:
                        break
                i = i-len(motif)
            else:
                i = i-1
    return motif_site,gaps2seq
#########################
######################################
def align2seqs(seq_a,seq_b,chr1,n,end,out_dir):#motif;gap;index;outdir
    fold = math.ceil(len(seq_b)/len(seq_a))
    if fold <=2:
        fold =2
    seq_a1 = seq_a * fold
    ########
    if out_dir[-1] == '/':
        out_dir = out_dir[:-1]
    dir1 = out_dir+'/gap_seq'+end
    dir2 = out_dir+'/gap_aln'+end
    #print(dir1)
    if os.path.exists(dir1):
        pass
    else:
        os.mkdir(dir1)
    if os.path.exists(dir2):
        pass
    else:
        os.mkdir(dir2)
    out_seq = dir1+'/'+chr1+'_seq'+str(n)
    with open(out_seq, "w+") as f1:
        f1.writelines(">A\n{0}\n>B\n{1}".format(seq_a1,seq_b))
    out_align = dir2+'/'+chr1+'_aln'+str(n)
    #os.system("chmod +777 "+ out_seq)
    #os.system("mafft --op 0.2 --quiet --clustalout "+out_seq+ " > "+ out_align)
    os.system("mafft --textmatrix matrix.txt --quiet --clustalout "+out_seq+ " > "+ out_align)
    #os.system("chmod +777 "+ out_align)
    my_file = open(out_align, "r")
    my_lines = my_file.readlines()
    seq_A = ''; seq_B= ''; symbol = ''
    for each1 in range(len(my_lines)):
        if my_lines[each1].startswith("A"):
            line1 = my_lines[each1][16:].strip("\n")
            seq_A += line1.upper()
        if my_lines[each1].startswith("B"):
            line2 = my_lines[each1][16:].strip("\n")
            seq_B += line2.upper()
            symbol += my_lines[each1+1][16:].strip('\n')
    #os.system('rm -r '+dir1)
    #os.system('rm -r '+dir2)
    return seq_A, seq_B, symbol
######################################

def BaseChange(symbol,motif,seq_A,seq_B):
    #print(sign)
    list1 = []
    i1 = 0
    while i1 < len(symbol):
        if '-' in symbol[i1:i1+len(motif)]:#insert
            x1 = 1
        else:
            x1 = 0
        symbols1 = symbol[i1:i1+len(motif)+x1]
        if len(symbols1) in [len(motif),len(motif)+1] and symbols1.count('*')  == len(motif)+x1-1:
            if list(set(symbol[i1+1:i1+len(motif)+x1+1])) !=1:
                symbols2 = list(set(symbols1))
                symbols2.remove('*')
                print(seq_B[i1:i1+len(motif)+x1])
                print(symbols1)
                sign = symbols2[0]
                i = symbols1.index(sign)
                #list1.append(i1+i)
                ###
                motif_tmp1 = seq_A[i1:i1+len(motif)+x1]
                #motif_tmp2 = motif_tmp1[ind:]+motif_tmp1[:ind]
                #print(motif_tmp1)
                #print(seq_B[i1:i1+len(motif)])
                print('!!',x1,sign,motif_tmp1,motif)
                if x1 ==0:
                    motifd = motif_tmp1*2#double
                    #print(motif_tmp1,motifd,motif)
                    ind = motifd.index(motif)
                else:
                    motif_tmp1 = list(motif_tmp1)
                    motif_tmp1.remove('-')
                    motif_tmp1 = ''.join(motif_tmp1)
                    motifd = motif_tmp1*2
                    ind = motifd.index(motif)
                #print(motif_tmp2)
                #print(x1,motif_tmp1,motifd,motif)
                if i >= ind:
                    i2 = i-ind+1
                    #print(i)
                    #print(i-ind+1)
                else:
                    i2 = len(motif)-x1-ind+i+1
                    #print(i)
                    #print(len(motif)-ind+i)
                ###
                i3 = i1+i
                symbols_tmp = symbol[:i3]
                for sym in symbol[:i3]:
                    if sym == ' ' or sym == '-':#del or ins
                        i3 =i3-1
                list1.append([i2,i3,seq_A[i1+i],seq_B[i1+i]])
                #rsite;abdis;b1;b2
                #在motif的相对位置；距离端粒末端的绝对距离
                i1 = i1+len(motif)+x1-1
        i1+=1
    print(seq_A)
    print(seq_B)
    print(list1)
    return list1

###################################
def symbol_count(seq_A,seq_B,motif,symbol):
    symbol = list(symbol)
    for i1 in range(len(symbol)):
        if seq_A[i1] == '-':#插入
            symbol[i1] = '-'
    for i1 in range(len(symbol)):
        if symbol[i1] == ' ':
            if seq_A[i1] !='-' and seq_B[i1] !='-':
                symbol[i1] = '.'#将颠换的space改为.
    #########
    all1 = BaseChange(symbol,motif,seq_A,seq_B)
    count_dc = {}
    for i1,i2,b1,b2 in all1:
        mut_type = '{0}>{1}'.format(b1,b2)
        if mut_type not in count_dc:
            count_dc[mut_type] = []
        count_dc[mut_type].append([i1,i2])
    #########
    return all1,count_dc

###################################
def IndelSeq(genomef,motif,end,out_dir):
    ####
    chr2n = {}#motif count
    chr2len = {}#telomere length assembled
    chr2gap = {}#chr:b1>b2,rsite,dis
    ####
    n1 = 1
    fas = SeqIO.parse(genomef,'fasta')
    tmp= []
    for k1 in fas:
        chr1 = k1.id
        if tmp == []:
            tmp.append(chr1[:2])
        if chr1[:2] == tmp[0] or len(chr1) <=5:
            seq1 = str(k1.seq).upper()
            motif_site,gaps2seq = gapfind(seq1,motif,end)
            motif_site = motif_site[:-1]
            chr2n[chr1] = len(motif_site)
            if end == 'left':
                lens = [a[1] for a in motif_site]
            if end == 'right':
                lens = [abs(a[0]) for a in motif_site]
            if len(motif_site) >=3 and max(lens) <= (len(motif)*2)*len(motif_site):
                chr2len[chr1] = max(lens)
            else:
                if len(motif_site)>=3:
                    chr2len[chr1] = len(motif)*len(motif_site)
                else:
                    chr2len[chr1] =0
            ###
            if chr2len[chr1] >0:
                for gaps in gaps2seq:
                    gseq = gaps2seq[gaps]
                    gs1,ge1 = [int(a) for a in gaps.split(',')]
                    count_dc = {}
                    if len(gseq) ==1:
                        if gseq == motif[0]:
                            rsite = 1
                        else:
                            rsite = len(motif)
                        count_dc['->'+gseq] = [[rsite,gs1+1]]
                    else: 
                        seq_A, seq_B, symbol=align2seqs(motif,gseq,chr1,n1,end,out_dir)
                        all1,count_dc=symbol_count(seq_A,seq_B,motif,symbol)
                        chr2n[chr1] = chr2n[chr1] + len(all1)
                        n1 +=1
                    if len(count_dc) >0:
                        if chr1 not in chr2gap:
                            chr2gap[chr1] = []
                        for mut in count_dc:
                            for rsite,dis in count_dc[mut]:
                                if len(gseq) ==1:
                                    chr2gap[chr1].append([mut,rsite,dis])
                                else:
                                    if end == 'left':
                                        dis1 = gs1+dis
                                    if end == 'right':
                                        dis1 = gs1+len(gseq)-dis
                                    #print(gs1,len(gseq),dis)
                                    #print('dis:',dis1)
                                    chr2gap[chr1].append([mut,rsite,dis1])
    return chr2n,chr2len,chr2gap

def OutputMut(info,motif,chr2n,chr2len,chr2gap,out_dir):
    ###h1
    h1 = open(os.path.join(out_dir,info+'_out1.txt'),'w')
    h1.write('\t'.join(['chr','motifn','SVn','TelLen','ab_dis'])+'\n')
    for chr1 in chr2n:
        motifn = chr2n[chr1]
        tellen = chr2len[chr1]
        if chr2len[chr1] >0:
            if chr1 in chr2gap:
                svn = len(chr2gap[chr1])
                dis_ls = [a[-1] for a in chr2gap[chr1]]
            else:
                svn = 0
                dis_ls = []
            dis_ls.sort()
            dis_ls = [str(a) for a in dis_ls]
            ###
            ###
            h1.write('\t'.join([chr1,str(motifn),str(svn),str(tellen),','.join(dis_ls)]))
            h1.write('\n')
    h1.close()
    ###h2
    count_dc1 = {}
    ls1 = []
    for chr1 in chr2gap:
        for mut,rsite,dis1 in chr2gap[chr1]:
            if rsite not in count_dc1:
                count_dc1[rsite] = []
                ls1.append(rsite-1)
            count_dc1[rsite].append(mut)
    h2 = open(os.path.join(out_dir,info+'_out2.txt'),'w')
    h2.write('motif:'+motif+'\n')
    ls1.sort()
    for rsite in ls1:
        h2.write(str(rsite+1)+'-'+motif[rsite]+':')
        h2.write(','.join(count_dc1[rsite+1])+'\n')
    h2.close()

def OutputMutfq(info,motif,chr2n,chr2len,chr2gap,out_dir):
    ###h1
    h1 = open(os.path.join(out_dir,info+'_out1.txt'),'w')
    h1.write('\t'.join(['motifn','SVn'])+'\n')
    motifn = sum([chr2n[a] for a in chr2n])
    svn = sum([len(chr2gap[a]) for a in chr2gap])
    h1.write('\t'.join([str(motifn),str(svn)]))
    h1.write('\n')
    h1.close()
    ###h2
    ###h2
    count_dc1 = {}
    ls1 = []
    for chr1 in chr2gap:
        for mut,rsite,dis1 in chr2gap[chr1]:
            if rsite not in count_dc1:
                count_dc1[rsite] = []
                ls1.append(rsite-1)
            count_dc1[rsite].append(mut)
    h2 = open(os.path.join(out_dir,info+'_out2.txt'),'w')
    h2.write('motif:'+motif+'\n')
    ls1.sort()
    for rsite in ls1:
        h2.write(str(rsite+1)+'-'+motif[rsite]+':')
        h2.write(','.join(count_dc1[rsite+1])+'\n')
    h2.close() 

