class HSP:
    def __init__(self, QNAME, QPOS, RNAME, RPOS):
        self.BIN_SIZE = 100
        self.SLIDING_SIZE = 100

        self.QNAME = QNAME
        self.QsPOS = QPOS
        self.QePOS = QPOS + self.BIN_SIZE - 1

        self.RNAME = RNAME
        self.RsPOS = RPOS
        self.RePOS = RPOS + self.BIN_SIZE - 1

        self.readN = 1

import math

def calc_idx(BIN_SIZE, SLIDING_SIZE, pos):
    maxIDX = math.floor(pos/SLIDING_SIZE)
    minIDX = math.ceil((pos - BIN_SIZE)/SLIDING_SIZE)

    return minIDX, maxIDX

BLOCK_W = 5000
BLOCK_H = 300
SLIDING = 100


from itertools import groupby
fin = open('query_100.sam')

fout = open('test.pos', 'w')

for line in fin:
    if line.startswith('@PG') == True:
        break

for QNAME, group1 in groupby(fin, lambda line: line.split('_')[0]):
    #print('----------------------------------------------------------', QNAME)
    count_DICT = {}
    hsp_LIST = []
    count = 0
    for key2, group2 in groupby(group1, lambda line: line.split('\t')[0]):
        #print('----------------------------------------------------------', key2)
        QNAME, QPOS = key2.split('_')
        QPOS = int(QPOS)
        for data in group2:
            FLAG, RNAME, RPOS, MAPQ, CIGAR = data.split('\t')[1:6]
            if RNAME == '*': continue
            RPOS = int(RPOS)
            FLAG = int(FLAG)
            STRAND = '+'
            if FLAG&16 == 16:
                STRAND = '-'
            hsp = [RNAME, STRAND, QPOS, RPOS, 0]
            hsp_LIST += [hsp]
            if not RNAME in count_DICT: count_DICT[RNAME] = {'P':{}, 'M':{}}

            intercept_p = -QPOS + RPOS
            intercept_m =  QPOS + RPOS

            #Slop: Plus ------------------------------------------------------------------------------------------
            intercept_p_minIDX = math.ceil((intercept_p - BLOCK_H)/SLIDING)
            intercept_p_maxIDX = math.floor((intercept_p)/SLIDING)

            intercept_m_minIDX = math.ceil((intercept_m - BLOCK_W)/SLIDING)
            intercept_m_maxIDX = math.floor((intercept_m)/SLIDING)

            for pIDX in range(intercept_p_minIDX, intercept_p_maxIDX + 1):
                for mIDX in range(intercept_m_minIDX, intercept_m_maxIDX + 1):
                    if not (pIDX, mIDX) in count_DICT[RNAME]['P']: 
                        count_DICT[RNAME]['P'][(pIDX, mIDX)]  = [hsp]
                    else:
                        count_DICT[RNAME]['P'][(pIDX, mIDX)] += [hsp]
            
            #Slop: Minus
            intercept_p_minIDX = math.ceil((intercept_p - BLOCK_W)/SLIDING)
            intercept_p_maxIDX = math.floor((intercept_p)/SLIDING)

            intercept_m_minIDX = math.ceil((intercept_m - BLOCK_H)/SLIDING)
            intercept_m_maxIDX = math.floor((intercept_m)/SLIDING)


            if intercept_p_minIDX > intercept_p_maxIDX:
                print(intercept_p_minIDX, intercept_p_maxIDX)
            for pIDX in range(intercept_p_minIDX, intercept_p_maxIDX + 1):
                for mIDX in range(intercept_m_minIDX, intercept_m_maxIDX + 1):
                    if not (pIDX, mIDX) in count_DICT[RNAME]['M']: 
                        count_DICT[RNAME]['M'][(pIDX, mIDX)]  = [hsp]
                    else:
                        count_DICT[RNAME]['M'][(pIDX, mIDX)] += [hsp]
    
    for RNAME, sub1Count_DICT in count_DICT.items():
        for slop, sub2Count_DICT in sub1Count_DICT.items():
            for (pIDX, mIDX), _hsp_LIST in sub2Count_DICT.items():
                if len(_hsp_LIST) >= BLOCK_W/SLIDING * 0.8:
                    for hsp in _hsp_LIST:
                        hsp[4] += 1

    count = 0
    for hsp in hsp_LIST:
        if hsp[4] > 0:
            count += 1
            fout.write('\t'.join(map(str, [QNAME] + hsp)) + '\n')
    fout.flush()
    print(QNAME, len(hsp_LIST), count)
fout.close()
