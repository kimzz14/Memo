import math

BLOCK_W = 1000*5 * 1.0
BLOCK_H = BLOCK_W * 0.05
SLIDING = BLOCK_H / 3

BLOCK_W = 3000
BLOCK_H = 100
SLIDING = 100

print('BLOCK_W:', BLOCK_W)
print('BLOCK_H:', BLOCK_H)
print('SLIDING:', SLIDING)

ReadSkip = 71

THRESHOLD = math.ceil(BLOCK_W / 2 / ReadSkip)

print(THRESHOLD) 

from itertools import groupby
fin = open('query_100.sam')

fout = open('test.pos', 'w')

for line in fin:
    if line.startswith('@PG') == True:
        break

for QNAME, group1 in groupby(fin, lambda line: line.split('_')[0]):
    #print('----------------------------------------------------------', QNAME)
    countP_DICT = {}
    countM_DICT = {}
    hsp_LIST = []
    count = 0
    for key2, group2 in groupby(group1, lambda line: line.split('\t')[0]):
        #print('----------------------------------------------------------', key2)
        QNAME, QPOS = key2.split('_')
        QPOS = int(QPOS)
        for data in group2:
            FLAG, RNAME, RPOS, MAPQ, CIGAR = data.split('\t')[1:6]
            if RNAME == '*': continue
            if not RNAME in countP_DICT: countP_DICT[RNAME] = {}
            if not RNAME in countM_DICT: countM_DICT[RNAME] = {}

            RPOS = int(RPOS)
            FLAG = int(FLAG)
            STRAND = '+'
            if FLAG&16 == 16:
                STRAND = '-'
            hsp = [RNAME, STRAND, QPOS, RPOS, 0]
            hsp_LIST += [hsp]
            
            intercept_p = -QPOS + RPOS
            intercept_m =  QPOS + RPOS 
            
            #Slop: Plus ------------------------------------------------------------------------------------------
            if STRAND == '+':
                intercept_p_minIDX = math.ceil((intercept_p - BLOCK_H)/SLIDING)
                intercept_p_maxIDX = math.floor((intercept_p)/SLIDING)

                intercept_m_minIDX = math.ceil((intercept_m - BLOCK_W)/SLIDING)
                intercept_m_maxIDX = math.floor((intercept_m)/SLIDING)

                for pIDX in range(intercept_p_minIDX, intercept_p_maxIDX + 1):
                    for mIDX in range(intercept_m_minIDX, intercept_m_maxIDX + 1):
                        if not (pIDX, mIDX) in countP_DICT[RNAME]: 
                            countP_DICT[RNAME][(pIDX, mIDX)]  = [hsp]
                        else:
                            countP_DICT[RNAME][(pIDX, mIDX)] += [hsp]
            
            #Slop: Minus------------------------------------------------------------------------------------------
            if STRAND == '-':
                intercept_p_minIDX = math.ceil((intercept_p - BLOCK_W)/SLIDING)
                intercept_p_maxIDX = math.floor((intercept_p)/SLIDING)

                intercept_m_minIDX = math.ceil((intercept_m - BLOCK_H)/SLIDING)
                intercept_m_maxIDX = math.floor((intercept_m)/SLIDING)

                for pIDX in range(intercept_p_minIDX, intercept_p_maxIDX + 1):
                    for mIDX in range(intercept_m_minIDX, intercept_m_maxIDX + 1):
                        if not (pIDX, mIDX) in countM_DICT[RNAME]: 
                            countM_DICT[RNAME][(pIDX, mIDX)]  = [hsp]
                        else:
                            countM_DICT[RNAME][(pIDX, mIDX)] += [hsp]

    print(QPOS)

    maxN = 0
    for RNAME, subCount_DICT in countP_DICT.items():
        for (pIDX, mIDX), _hsp_LIST in subCount_DICT.items():
            if len(_hsp_LIST) >= THRESHOLD * 0.8:
                maxN = max(maxN, len(_hsp_LIST))
                for hsp in _hsp_LIST:
                    hsp[4] += 1
    
    for RNAME, subCount_DICT in countM_DICT.items():
        for (pIDX, mIDX), _hsp_LIST in subCount_DICT.items():
            if len(_hsp_LIST) >= THRESHOLD * 0.8:
                maxN = max(maxN, len(_hsp_LIST))
                for hsp in _hsp_LIST:
                    hsp[4] += 1

    count = 0
    test_DICT = {}
    for hsp in hsp_LIST:
        if hsp[4] > 0:
            count += 1
            fout.write('\t'.join(map(str, [QNAME] + hsp)) + '\n')
            if not hsp[0] in test_DICT: test_DICT[hsp[0]] = [1000000000, 0, 1000000000, 0]
            test_DICT[hsp[0]][0] = min(test_DICT[hsp[0]][0], hsp[2])
            test_DICT[hsp[0]][1] = max(test_DICT[hsp[0]][1], hsp[2])

            test_DICT[hsp[0]][2] = min(test_DICT[hsp[0]][2], hsp[3])
            test_DICT[hsp[0]][3] = max(test_DICT[hsp[0]][3], hsp[3])
    fout.flush()
    print(QNAME, len(hsp_LIST), count, maxN)

    for RNAME ,count in test_DICT.items():
        print(QNAME, RNAME, count, count[1] - count[0], count[3] - count[2])
    break
fout.close()
