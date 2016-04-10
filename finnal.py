

def orf_finder(seq, frame):
    i1 = []
    i2 = []
    seq = seq[frame-1:]
    counts = len(seq) // 3
    for c in range(counts):
        index = c*3
        tmp = seq[index: index+3]
        if tmp == 'ATG':
            i1.append(index)
        elif tmp in ['TGA', 'TAG', 'TAA']:
            i2.append(index)

    result = []
    for i in i1:
        for j in i2:
            if j > i:
                p1 = i
                p2 = j+3
                tmp = {
                    'len': p2-p1,
                    'p1': p1,
                    'p2': p2,
                    'orf': seq[p1:p2]
                }
                result.append(tmp)
                break
    return result

def concat_list(orf_info):
    lists = []
    for i in orf_info:
        lists += i
    return lists

def seq_counter(seq, length, result):
    total_len = len(seq)
    for n in range(length, total_len):
        tmp = seq[n-length: n]
        if tmp not in result:
            result[tmp] = 1
        else:
            result[tmp] += 1


seqs = {}
with open('dna2.fasta') as f:
    for line in f:
        line = line.rstrip()
        if line[0] == '>':
            words = line.split()
            name = words[0][1:]
            seqs[name] = ''
        else:
            seqs[name] = seqs[name] + line

# find ORF
seq_values = seqs.values()
orf1 = {name: orf_finder(s, 1) for name, s in seqs.items()}
orf2 = {name: orf_finder(s, 2) for name, s in seqs.items()}
orf3 = {name: orf_finder(s, 3) for name, s in seqs.items()}

dict_list = []

# dict list 
list1 = concat_list(orf1.values())
list2 = concat_list(orf2.values())
list3 = concat_list(orf3.values())
# orf lens
orf1_lens = [v['len'] for v in list1]
orf2_lens = [v['len'] for v in list2]
orf3_lens = [v['len'] for v in list3]

# question 1
print 'Q1:', len(seqs)
# question 2
seq_lens = [len(s) for s in seq_values]
print 'Q2:', max(seq_lens)
# question 3
print 'Q3:', min(seq_lens)
# question 4
print 'Q4:', max(orf2_lens)
# question 5
postions_3 = [v['p1'] for v in list3]
logest_index = orf3_lens.index(max(orf3_lens))
# plus 3 for the first postion index is 0 in programing
print 'Q5:', postions_3[logest_index] + 3
# question 6
total_orf_lens = orf1_lens + orf2_lens + orf3_lens
print 'Q6:', max(total_orf_lens)
# question 7
name = 'gi|142022655|gb|EQ086233.1|16'
s_orfs = orf1[name] + orf2[name] + orf3[name]
s_lens = [v['len'] for v in s_orfs]
print 'Q7:', max(s_lens)
# question 8
result8 = {}
for s in seq_values:
    seq_counter(s, 6, result8)
print 'Q8:', max(result8.values())
# question 9
result9 = {}
for s in seq_values:
    seq_counter(s, 12, result9)
max_num = max(result9.values())
print 'Q9:', result9.values().count(max_num)
# question 10
result10 = {}
for s in seq_values:
    seq_counter(s, 7, result10)
max_num = 0
max_seq = ''
for s, num in result10.items():
    if num > max_num:
        max_num = num
        max_seq = s
print 'Q10:', max_seq