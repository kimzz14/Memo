def TmCalculator(Primer): #TmCalculator v1.0.0, compiled 2016/04/11 12:38:00
    import math
    Primer = Primer.upper()
    bh = ""
    sumH = 0.0
    sumS = 0.0
    concSalt = 50.0
    concPrimer = 200.0
    h = {}
    s = {}
    h['A'] = {}
    h['C'] = {}
    h['G'] = {}
    h['T'] = {}
    s['A'] = {}
    s['C'] = {}
    s['G'] = {}
    s['T'] = {}
    h['A']['A'] = 9.1
    h['A']['C'] = 6.5
    h['A']['G'] = 7.8
    h['A']['T'] = 8.6
    h['C']['A'] = 5.8
    h['C']['C'] = 11
    h['C']['G'] = 11.9
    h['C']['T'] = 7.8
    h['G']['A'] = 5.6
    h['G']['C'] = 11.1
    h['G']['G'] = 11
    h['G']['T'] = 6.5
    h['T']['A'] = 6
    h['T']['C'] = 5.6
    h['T']['G'] = 5.8
    h['T']['T'] = 9.1
    s['A']['A'] = 24
    s['A']['C'] = 17.3
    s['A']['G'] = 20.8
    s['A']['T'] = 23.9
    s['C']['A'] = 12.9
    s['C']['C'] = 26.6
    s['C']['G'] = 27.8
    s['C']['T'] = 20.8
    s['G']['A'] = 13.5
    s['G']['C'] = 26.7
    s['G']['G'] = 26.6
    s['G']['T'] = 17.3
    s['T']['A'] = 16.8
    s['T']['C'] = 13.5
    s['T']['G'] = 12.9
    s['T']['T'] = 24
    isFirst = True
    for ch in Primer:
        if isFirst:
            isFirst = False
            bh = ch
            continue
        sumH += h[bh][ch]
        sumS += s[bh][ch]
        bh = ch

    if Primer[0] == 'G' or Primer[0] == 'C':
        sumH += 0.25/2
        sumS += -0.62/2
    else:
        sumH += 0.7/2
        sumS += 0.62/2
    if Primer[-1] == 'G' or Primer[-1] == 'C':
        sumH += 0.25/2
        sumS += -0.62/2
    else:
        sumH += 0.7/2
        sumS += 0.62/2
        bh = ch

    return -1000*sumH/(-1*(sumS+16.8)+1.987*math.log(concPrimer/4000000000))-273.15+16.6*math.log(concSalt/1000)/math.log(10);


fin = open('PrimerLIST.txt')
fout = open('PrimerLIST.Tm.txt', 'w')

fout.write('\t'.join(map(str,['fPrimer','rPrimer','fLen','rLen','fTm','rTm','fGC','rGC'])) + '\n')
for line in fin:
    fPrimer, rPrimer = line.rstrip('\n').upper().split('\t')
    fTm  = TmCalculator(fPrimer)
    rTm  = TmCalculator(rPrimer)
    fLen = len(fPrimer)
    rLen = len(rPrimer)
    fGC = float(fPrimer.count('G') + fPrimer.count('C')) / fLen
    rGC = float(rPrimer.count('G') + rPrimer.count('C')) / rLen
    fout.write('\t'.join(map(str,[fPrimer,rPrimer,fLen,rLen,fTm,rTm,fGC,rGC])) + '\n')
fout.close()
