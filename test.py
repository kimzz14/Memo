class HSP:
    def __init__(self,queryPos, sbjctPos):
        self.sQueryPos = queryPos
        self.sSbjctPos = sbjctPos
        self.eQueryPos = queryPos
        self.eSbjctPos = sbjctPos
    
    def isNextQuery(self,queryPos):
        if self.eQueryPos + 1 != queryPos: 
            return False
        return True
    def isNextSbjct(self,sbjctPos):
        if self.eSbjctPos + 1 != sbjctPos: 
            return False
        return True

    def add(self,queryPos, sbjctPos):
        self.eQueryPos = queryPos
        self.eSbjctPos = sbjctPos

    def text(self):
        return '\t'.join(map(str,[self.eQueryPos - self.sQueryPos + 1, self.sQueryPos, self.eQueryPos, self.sSbjctPos, self.eSbjctPos]))
        #return '[' + str(self.eQueryPos - self.sQueryPos + 1) + ']' + '\t' + str(self.sQueryPos) + ' ~ ' + str(self.eQueryPos) + '   ' + str(self.sSbjctPos) + ' ~ ' + str(self.eSbjctPos)


HSP_DICT = {}

fin = open('all.sam')

for line in fin:
    if line.startswith('@PG') == True:
        break

for lineIDX, line in enumerate(fin):
    data_LIST = line.rstrip('\n').split('\t')
    queryID = data_LIST[0]

    query, pos_LIST = queryID.split('_')
    queryPos = int(pos_LIST.split('-')[0])


    if not query in HSP_DICT: HSP_DICT[query] = {}


    match_LIST = []

    match = []
    sbjct = data_LIST[2]
    if sbjct == '*': continue

    pos = int(data_LIST[3])
    flag = int(data_LIST[1])
    if flag&16 == 16:
        pos *= -1

    match_LIST += [(sbjct, pos)]

    for data in data_LIST:
        if data.startswith('XA') == True:
            for match in data[5:].split(';')[:-1]:
                sbjct, pos, cigar, tmp = match.split(',')
                pos = int(pos)
                match_LIST += [(sbjct, pos)]

    for match in match_LIST:
        sbjct, sbjctPos = match
        if not sbjct in HSP_DICT[query]: 
            HSP_DICT[query][sbjct] = [HSP(queryPos, sbjctPos)]
        else:
            hsp_LIST = []
            for hsp in HSP_DICT[query][sbjct]:
                if hsp.isNextQuery(queryPos) == False:
                    if hsp.eQueryPos - hsp.sQueryPos + 1 > 10:
                        print(query + '\t' +  sbjct + '\t' +  hsp.text())
                else:
                    hsp_LIST += [hsp]
                    if hsp.isNextSbjct(sbjctPos) == True:
                        hsp.add(queryPos, sbjctPos)
            if len(hsp_LIST) == 0:
                hsp_LIST += [HSP(queryPos, sbjctPos)]
            HSP_DICT[query][sbjct] = hsp_LIST
            
for query, item in HSP_DICT.items():
    for sbjct, hsp_LIST in item.items():
        for hsp in hsp_LIST:
            print(query + '\t' +  sbjct + '\t' +  hsp.text())
