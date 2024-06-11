from itertools import groupby
import gzip

target_LIST = []
target_DICT = {}


fin = open('IsolatedSNP.Exon.61mer_LIST')
for line in fin:
    record = line.rstrip('\n').split('\t')
    seqName, pos, prev, next = record
    pos = int(pos)
    
    if not seqName in target_DICT: target_DICT[seqName] = {}
    
    target_LIST += [record]
    target_DICT[seqName][pos] = record



fin = gzip.open('pooled.RGsorted.dedupped.HaplotypeCaller.all.vcf.gz', 'rt')
fout = open('IsolatedSNP_LIST', 'w')
for line in fin:
    if line.startswith('#CHROM') == True:
        legend_LIST = line.rstrip('\n').split('\t')
        break

genoIDX = legend_LIST.index('FORMAT') + 1
print(genoIDX)
countN = 0
for line in fin:
    data_LIST = line.rstrip('\n').split('\t')
    seqName = data_LIST[0]
    pos = int(data_LIST[1])

    if not seqName in target_DICT: continue
    if not pos in target_DICT[seqName]: continue
    

    missCallN = 0
    totalCallN = 0
    for data in data_LIST[genoIDX:]:
        totalCallN += 1
        genotype = data.split(':')[0]

        if genotype == './.':
            missCallN += 1

    callRate = 1 - float(missCallN)/totalCallN
    #print(missCallN, totalCallN, callRate)

    record = target_DICT[seqName][pos]
    record += [missCallN, totalCallN, callRate]

    if callRate > 0.95:
        countN += 1
print(countN)

fout = open('IsolatedSNP.Exon.61mer.callRate_LIST', 'w')

for target in target_LIST:
    fout.write('\t'.join(map(str, target)) + '\n')
fout.close()
