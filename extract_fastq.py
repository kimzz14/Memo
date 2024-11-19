readID_DICT = {}
contigID_DICT = {}
fin = open('hifiasm.asm.bp.p_ctg.gfa')

for line in fin:
    if line[0] != 'A': continue

    data_LIST = line.rstrip('\n').split('\t')
    #print(data_LIST)
    contigID = data_LIST[1]
    readID = data_LIST[4]

    if readID in readID_DICT and readID_DICT[readID] != contigID: 
        print('bug:', readID, contigID, readID_DICT[readID])

    readID_DICT[readID] = contigID

    if not contigID in contigID_DICT: contigID_DICT[contigID] = 0
fin.close()


contigID_LIST = sorted(list(contigID_DICT.keys()))

import gzip

fin = gzip.open('SRR19088064_subreads.fastq.gz', 'rt')

fout_DICT = {}
for contigID in contigID_LIST + ['ETC']:
    fout_DICT[contigID] = gzip.open('result/' + contigID + '.fastq.gz', 'wt')

for lineIDX, line in enumerate(fin):
    if lineIDX%4 == 0:
        readID = line.rstrip('\n').split(' ')[0].split('\t')[0][1:]

        if readID in readID_DICT:
            contigID = readID_DICT[readID]
        else:
            contigID = 'ETC'
        #print(readID, contigID)


    fout_DICT[contigID].write(line)

for contigID in contigID_LIST + ['ETC']:
    fout_DICT[contigID].close()
