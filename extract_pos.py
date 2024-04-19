class HSP:
    def __init__(self, QNAME, QPOS, RNAME, RPOS):
        self.BIN_SIZE = 100
        self.SLIDING_SIZE = 71

        self.QNAME = QNAME
        self.QsPOS = QPOS
        self.QePOS = QPOS + self.BIN_SIZE - 1

        self.RNAME = RNAME
        self.RsPOS = RPOS
        self.RePOS = RPOS + self.BIN_SIZE - 1

    def Q_isNext(self, QNAME, QPOS):
        if self.QNAME == QNAME:
            if self.QePOS + 1 == QPOS + self.BIN_SIZE - self.SLIDING_SIZE:
                return True
        return False

    def R_isNext(self, RNAME, RPOS):
        if self.RNAME == RNAME:
            if self.RePOS + 1 == RPOS + self.BIN_SIZE - self.SLIDING_SIZE:
                return True
        return False

    def add(self, QNAME, QPOS, RNAME, RPOS):
        self.QePOS = QPOS + self.BIN_SIZE - 1
        self.RePOS = RPOS + self.BIN_SIZE - 1

    def info(self):
        return [self.QNAME, self.RNAME, self.QePOS - self.QsPOS + 1, self.RePOS - self.RsPOS + 1, self.QsPOS, self.QePOS, self.RsPOS, self.RePOS]

from itertools import groupby
fin = open('head.sam')

for line in fin:
    if line.startswith('@PG') == True:
        break

nHSP_LIST = []

for key, group in groupby(fin, lambda line: line.split('\t')[0]):
    #print('----------------------------------------------------------', key)
    QNAME, QPOS = key.split('_')
    QPOS = int(QPOS)

    cHSP_LIST = nHSP_LIST
    nHSP_LIST = []
    dHSP_LIST = []
    for cHSP in cHSP_LIST:
        if cHSP.Q_isNext(QNAME, QPOS) == True:
            nHSP_LIST += [cHSP]
        else:
            dHSP_LIST += [cHSP]

    for hsp in dHSP_LIST:
        #print("aaa", hsp.RNAME, hsp.QePOS - hsp.QsPOS + 1, hsp.RePOS - hsp.RsPOS + 1, hsp.QsPOS, hsp.QePOS, hsp.RsPOS, hsp.RePOS)
        print('\t'.join(map(str, hsp.info())))
    
    dHSP_LIST = nHSP_LIST
    nHSP_LIST = []
    for data in group:
        #print('[group] dHSP_LIST:', len(dHSP_LIST))
        #print('[group] nHSP_LIST:', len(nHSP_LIST))
        cHSP_LIST = dHSP_LIST
        dHSP_LIST = []

        FLAG, RNAME, RPOS, MAPQ, CIGAR = data.split('\t')[1:6]

        FLAG = int(FLAG)
        if FLAG&16 == 16:
            RPOS = int('-' + RPOS)
        else:
            RPOS = int(RPOS)

        #print(RNAME, RPOS, CIGAR)
        
        isAdd = False
        for cHSP in cHSP_LIST:
            if cHSP.R_isNext(RNAME, RPOS) == True:
                cHSP.add(QNAME, QPOS, RNAME, RPOS)
                nHSP_LIST += [cHSP]
                isAdd = True
            else:
                dHSP_LIST += [cHSP]
        if isAdd == False:
            nHSP_LIST += [HSP(QNAME, QPOS, RNAME, RPOS)]

    #dHSP_LIST
    #print('[out] dHSP_LIST:', len(dHSP_LIST))
    #print('[out] nHSP_LIST:', len(nHSP_LIST))

    for hsp in dHSP_LIST:
        #print("bbb", hsp.info())
        print('\t'.join(map(str, hsp.info())))


for hsp in nHSP_LIST:
    #print("ccc", nHSP.info())
    print('\t'.join(map(str, hsp.info())))
