
import os
import math

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

######################################
def align2seqs(seq_a,seq_b,n,out_dir):
    fold = math.ceil(len(seq_b)/len(seq_a))
    seq_a1 = seq_a * fold
    dir1 = out_dir+'/tmp'
    if os.path.exists(dir1):
        pass
    else:
        os.mkdir(dir1)
    out_seq = out_dir + '/tmp/temp_seq'+str(n)
    with open(out_seq, "w+") as f1:
        f1.writelines(">A\n{0}\n>B\n{1}".format(seq_a1,seq_b))
    out_align = out_dir + '/tmp/temp_align'+str(n)
    os.system("chmod +777 "+ out_seq)
    os.system("sudo mafft --quiet --clustalout --auto "+out_seq+ " > "+ out_align)
    os.system("chmod +777 "+ out_align)
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
    return seq_A, seq_B, symbol

###################################
def BaseChange(symbol,sign):
    list1 = []
    if sign == ' ' or sign == '-' or sign == '.':
        for i1 in range(len(symbol)):
            if symbol[i1] == sign:  # 缺失：sign=' '；插入：sign = '-'（单碱基）
                if i1 == 0:
                    if symbol[i1 + 1] != sign:
                        list1.append(i1)
                if i1 == len(symbol) - 1:
                    if symbol[i1 - 1] != sign:
                        list1.append(i1)
                if i1 > 0 and i1 < len(symbol) - 1:
                    if symbol[i1 - 1] != sign and symbol[i1 +1] != sign:
                        list1.append(i1)
    #if sign == '.':#突变(多碱基)
        #for i1 in range(len(symbol)):
            #if symbol[i1] == '.':
                #list1.append(i1)
    if sign == '-':#记录插入的位置（由于多碱基插入会被忽略）
        list0 = []
        for i1 in range(len(symbol)):
            if symbol[i1] == '-':
                list0.append(i1)
        list1.append(list0)
    return list1
###################################
def symbol_count(seq_A,seq_B,motif1,symbol,index):
    symbol = list(symbol)
    for i1 in range(len(symbol)):
        if seq_A[i1] == '-':#插入
            symbol[i1] = '-'#将插入的space改为-
        if symbol[i1] == ' ':
            if seq_A[i1] !='-' and seq_B[i1] !='-':
                symbol[i1] = '.'#将颠换的space改为.
    symbol = ''.join(symbol)
    ist_ls = BaseChange(symbol, '-')  # 插入位置
    mut_ls = BaseChange(symbol, '.')  # 突变位置
    del_ls = BaseChange(symbol, ' ')  # 缺失位置
    all = ist_ls[:-1] + mut_ls + del_ls
    all.sort()
    ind = 0
    count_dc = {}
    for i1 in range(len(symbol)):
        if i1 in ist_ls[-1]:
            ind += 1
        if i1 in all:
            site = (i1 + 1 - ind) % len(motif1)
            if site == 0:
                site = len(motif1)
            # print(site+index)
            std_site = (index + site) % len(motif1)
            if std_site == 0:
                std_site = len(motif1)
            mut_type = "{0}>{1}".format(seq_A[i1], seq_B[i1])
            # print(str(i1)+':'+str(site)+':'+str(std_site)+':'+mut_type)
            if mut_type not in count_dc:
                count_dc[mut_type] = [std_site]
            else:
                count_dc[mut_type].append(std_site)
    return count_dc

###################################
def IndelSeq(input_file,motif):
    chr_motif = {}
    chr_seq = {}
    with open(input_file) as f1:
        for lines in f1.readlines()[1:]:
            line = lines.strip().split('\t')
            kmer1 = line[2]
            out1 = compare_kmer(motif, kmer1)
            if out1 == 'same':
                if line[4] == '0':
                    pass
                else:
                    chr_motif[line[0]] = line[2]
                    ls1 = []
                    for each1 in line[6].split(';')[:-1]:
                        ls1.append(each1)
                    chr_seq[line[0]] = ls1
    return chr_motif, chr_seq
###################################
def OutputMut(chr_motif, chr_seq, motif, info, out_dir):
    n1 = 1
    h1 = open(os.path.join(out_dir,info+'_out1'), 'w+')
    all_counts = {}
    for k1 in chr_motif:
        my_motif = chr_motif[k1]
        motif1 = motif * 2  # identify index to standard motif
        index1 = motif1.find(my_motif)
        seq_ls = chr_seq[k1]
        h1.write(k1 + '\t')  # chr
        # print(k1)
        for my_seq in seq_ls:
            if len(my_seq) == 1:
                count_dict1 = {}
                mut_type = '->' + my_seq
                site0 = (index1 + len(my_motif)) % len(my_motif)
                count_dict1[mut_type] = [site0]
            else:
                seq_A, seq_B, symbols = align2seqs(my_motif, my_seq, n1, out_dir)
                n1 += 1
                count_dict1 = symbol_count(seq_A, seq_B, my_motif,symbols, index1)
            if len(count_dict1) != 0:
                for k2 in count_dict1:
                    h1.write(k2 + ':')
                    if k2 not in all_counts:
                        all_counts[k2] = []
                    for j0 in count_dict1[k2]:
                        all_counts[k2].append(j0)
                    if len(count_dict1[k2]) == 1:
                        h1.write(str(count_dict1[k2][0]) + ';')
                    else:
                        for k3 in range(len(count_dict1[k2]) - 1):
                            h1.write(str(count_dict1[k2][k3]) + ',')
                        h1.write(';')
        h1.write('\n')
    h1.close()
    ###############
    h2 = open(os.path.join(out_dir, info+'_out2'), 'w+')
    h2.write('type:site:'+motif+'\n')
    for k1 in all_counts:
        h2.write(k1 + ':')
        for k2 in all_counts[k1]:
            h2.write(str(k2) + ',')
        h2.write('\n')
    h2.close()
    ###############
    all_count1 = {}
    for k1 in all_counts:
        for k2 in all_counts[k1]:
            if k2 not in all_count1:
                all_count1[k2] = [k1]
            else:
                all_count1[k2].append(k1)
    h3 = open(os.path.join(out_dir,info+'_out3'), 'w+')
    h3.write('site:type:'+motif+'\n')
    for k1 in all_count1:
        h3.write(str(k1)+'-'+motif[k1-1] + ':')
        for k2 in all_count1[k1]:
            h3.write(k2 + ',')
        h3.write('\n')
    h3.close()
    ################
    return 0
###################################



