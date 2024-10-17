inputDir = 'input'

ref_DICT = {}
seqName_LIST = []
fin = open('iwgsc_refseq_v2.1/ref.fa')
for lineIDX, line in enumerate(fin):
    if lineIDX > 5000: break
    if line.startswith('>') == True:
        seqName = line[1:].rstrip('\n').split(' ')[0].split('\t')[0]
        seqName_LIST += [seqName]
        ref_DICT[seqName] = []
    else:
        sequence = line.rstrip('\n')
        ref_DICT[seqName] += [sequence]

    
fin.close()

for seqName in seqName_LIST:
    ref_DICT[seqName] = ''.join(ref_DICT[seqName])


fout_command = open('command.sh', 'w')

fin = open('iwgsc_refseq_v2.1/ori/annotation/tmp/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3')
fin.readline()
for line in fin:
    seqName, source, feature, start, end, score, strand, frame, description = line.rstrip('\n').split('\t')
    if feature != 'gene': continue

    geneID = description.split(';')[0].split('=')[1]

    #print(geneID, start, end)

    fout = open(inputDir + '/' + geneID + '.fa', 'w')
    fout.write('>' + geneID + ':' + start + '-' + end + '\n')
    fout.write(ref_DICT[seqName][int(start) - 1: int(end)])
    fout.close()
    fout_command.write('sh blastn ' + geneID + '\n')
    break
fout_command.close()
