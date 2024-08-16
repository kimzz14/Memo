import gzip, math


splitSize = 2000

fin = gzip.open('../Keumgang.03cell.hifi.fastq.gz', 'rt')

fout = open('Keumgang.03cell.hifi.split2k.fastq', 'w')

readN = 0
for lineIDX, line in enumerate(fin):
    if lineIDX%4 == 0:
        readN += 1
        readID = 'Keumgang-' + str(readN)
    elif lineIDX%4 == 1:
        sequence = line.rstrip('\n')
        seqLen = len(sequence)
    elif lineIDX%4 == 2:
        pass

    elif lineIDX%4 == 3:
        quality = line.rstrip('\n')

        seqLen = len(sequence)
        if seqLen < 8000: continue

        for subIDX in range(0, math.ceil(seqLen/splitSize)):
            fout.write('@' + readID + '-' + str(seqLen) + '-' + str(splitSize*subIDX) + '\n')
            fout.write(sequence[splitSize*subIDX: splitSize*(subIDX + 1)] + '\n')
            fout.write('+' + '\n')
            fout.write(quality[splitSize*subIDX: splitSize*(subIDX + 1)] + '\n')
    else:
        print('bug')
