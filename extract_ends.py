revNucl_DICT = {}
revNucl_DICT['A'] = 'T'
revNucl_DICT['T'] = 'A'
revNucl_DICT['G'] = 'C'
revNucl_DICT['C'] = 'G'

def reverseComplementary(sequence):
    result = []
    for nucl in sequence[::-1]:
        result += revNucl_DICT[nucl]
    return ''.join(result)



flankSize = 30000

ref_DICT = {}
seqName_LIST = []

fin = open('hifiasm.asm.hic.p_ctg.fa')

for line in fin:
    if line.startswith('>') == True:
        seqName = line.rstrip('\n')[1:]
        ref_DICT[seqName] = []
        seqName_LIST += [seqName]
    else:
        sequence = line.rstrip('\n')
        ref_DICT[seqName] += [sequence]
fin.close()

for seqName in seqName_LIST:
    ref_DICT[seqName] = ''.join(ref_DICT[seqName])

fout_command = open('run01_blastn.sh', 'w')
for seqName in seqName_LIST:
    sequence = ref_DICT[seqName]

    #forward
    sPos = 0
    ePos = min(flankSize, len(sequence))
    fseqName = seqName + '-' + 'f' + '-' + str(sPos + 1) + '-' + str(ePos)
    fsequence = sequence[sPos:ePos]
    fout = open('input/' + fseqName + '.fa', 'w')
    fout.write('>' + fseqName + '\n')
    fout.write(fsequence + '\n')
    fout.close()
    fout_command.write(f'blastn -word_size 100 -perc_identity 98 -soft_masking false -dust no -query input/{fseqName}.fa -db hifiasm.asm.hic.p_ctg.blastDB/ref.fa -outfmt 6 > result/{fseqName}.blastn_m8\n')

    #reverse
    ePos = len(sequence)
    sPos = max(ePos - flankSize, 1)
    rseqName = seqName + '-' + 'r' + '-' + str(sPos +1) + '-' + str(ePos)
    rsequence = sequence[sPos:ePos]
    fout = open('input/' + rseqName + '.fa', 'w')
    fout.write('>' + rseqName + '\n')
    fout.write(rsequence + '\n')
    fout.close()
    fout_command.write(f'blastn -word_size 100 -perc_identity 98 -soft_masking false -dust no -query input/{rseqName}.fa -db hifiasm.asm.hic.p_ctg.blastDB/ref.fa -outfmt 6 > result/{rseqName}.blastn_m8\n')

fout_command.close()
