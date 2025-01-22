import glob

ref_DICT = {}
fin = open('hifiasm.asm.hic.p_ctg.fa.fai')
for line in fin:
    seqid, seqlen = line.rstrip('\n').split('\t')[0:2]
    ref_DICT[seqid] = int(seqlen)
fin.close()

for fileName in sorted(glob.glob('result/*_m8')):
    fin = open(fileName)
    for line in fin:
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.rstrip('\n').split('\t')
        if qseqid.split('-')[0] == sseqid: continue
        sstart = int(sstart)
        send = int(send)
        qseqid, qstrand, _qstart, _qend = qseqid.split('-')

        qstart = int(qstart) + int(_qstart) - 1
        qend = int(qend) + int(_qstart) - 1


        if sstart < send:
            strand = '+'
        else:
            strand = '-'

        length = int(length)
        pident = float(pident)

        
        sLen = ref_DICT[sseqid]

        srate = sstart / sLen
        
        if pident > 99.9 and length > 10000:
            if qstrand == 'r' and strand == '+' and srate < 0.2:
                print(qseqid, sseqid, qstrand, strand, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sLen, f'{srate:0.5f}')
            if qstrand == 'r' and strand == '-' and srate > 0.8:
                print(qseqid, sseqid, qstrand, strand, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sLen, f'{srate:0.5f}')

            if qstrand == 'f' and strand == '+' and srate > 0.8:
                print(qseqid, sseqid, qstrand, strand, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sLen, f'{srate:0.5f}')
            
            if qstrand == 'f' and strand == '-' and srate < 0.2:
                print(qseqid, sseqid, qstrand, strand, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, sLen, f'{srate:0.5f}')
