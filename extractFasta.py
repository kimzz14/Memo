class GENO:
    def __init__(self):
        self.sampleID_LIST = []
        self.geno_DICT = {}
        self.probeID_DICT = {}
        self.probeID_LIST = []

        self.fileN = 0

    def read_geno(self, prefix):
        self.fileN += 1
        fin = open('/pbi-acc2/mwKim/04_GWAS/02.GAPIT3/DataSet/genotype/' + prefix + '/' + 'genotype_calls.standard.txt')

        legend_LIST = fin.readline().rstrip('\n').split('\t')

        sampleID_LIST = legend_LIST[1:]
        probeID_DICT = {}

        for sampleID in sampleID_LIST:
            if not sampleID in self.geno_DICT:
                self.geno_DICT[sampleID] = {}
            else:
                print("duplicated sampleID:", sampleID)
        count = 0
        for line in fin:
            count += 1
            if count > 10000: break
            data_LIST = line.rstrip('\n').split('\t')
            fields = {key: value for key, value in zip(legend_LIST, data_LIST)}
            probeID = fields['probeset_id']

            if not probeID in probeID_DICT:
                probeID_DICT[probeID]  = 1
            else:
                probeID_DICT[probeID] += 1
                print("duplicated probeID:", probeID)

            for sampleID in sampleID_LIST:
                self.geno_DICT[sampleID][probeID] = fields[sampleID]
        fin.close()

        self.sampleID_LIST += sampleID_LIST
        for probeID, probeN in probeID_DICT.items():
            if not probeID in self.probeID_DICT:
                self.probeID_DICT[probeID]  = probeN
            else:
                self.probeID_DICT[probeID] += probeN

        print('Read fileName: {0}, probeN: {1}, sampleN: {2}'.format(prefix, len(probeID_DICT.keys()), len(sampleID_LIST)))

    def make_fasta(self, prefix, probeID_LIST, call_rate_threshold=0, sampleID_LIST=None):
        if sampleID_LIST is None:
            sampleID_LIST = self.sampleID_LIST

        passProbeID_DICT = {}
        for probeID in probeID_LIST:
            if not probeID in self.probeID_DICT: continue
            if self.probeID_DICT[probeID] != self.fileN: continue

            noCallN = 0
            for sampleID in sampleID_LIST:
                if not probeID in self.geno_DICT[sampleID]:
                    noCallN += 1
                elif self.geno_DICT[sampleID][probeID] == 'NN':
                    noCallN += 1

            call_rate = 1 - (noCallN / len(sampleID_LIST))

            if call_rate >= call_rate_threshold:
                passProbeID_DICT[probeID] = call_rate

        fout = open('{0}_S{1:05}_P{2:05}_CR{3:02.0f}.probe'.format(prefix, len(sampleID_LIST), len(passProbeID_DICT.keys()), call_rate_threshold * 100), 'w')
        for probeID in probeID_LIST:
            if probeID in passProbeID_DICT:
                fout.write(probeID + '\n')
        fout.close()

        fout = open('{0}_S{1:05}_P{2:05}_CR{3:02.0f}.fa'.format(prefix, len(sampleID_LIST), len(passProbeID_DICT.keys()), call_rate_threshold * 100), 'w')
        for sampleID in sampleID_LIST:
            fout.write('>' + sampleID + '\n')
            geno_LIST = []
            for probeID in probeID_LIST:
                if not probeID in passProbeID_DICT: continue
                geno_LIST += self.geno_DICT[sampleID][probeID]

            fout.write(''.join(geno_LIST) + '\n')
        fout.close()

        print('Write fileName: {0}, probeN: {1}, sampleN: {2}'.format(prefix, len(passProbeID_DICT.keys()), len(sampleID_LIST)))

    def make_geno(self, prefix, probeID_LIST, call_rate_threshold=0, sampleID_LIST=None):
        if sampleID_LIST is None:
            sampleID_LIST = self.sampleID_LIST

        passProbeID_DICT = {}
        for probeID in probeID_LIST:
            if not probeID in self.probeID_DICT: continue
            if self.probeID_DICT[probeID] != self.fileN: continue

            noCallN = 0
            for sampleID in sampleID_LIST:
                if not probeID in self.geno_DICT[sampleID]:
                    noCallN += 1
                elif self.geno_DICT[sampleID][probeID] == 'NN':
                    noCallN += 1

            call_rate = 1 - (noCallN / len(sampleID_LIST))

            if call_rate >= call_rate_threshold:
                passProbeID_DICT[probeID] = call_rate

        fout = open('{0}_S{1:05}_P{2:05}_CR{3:02.0f}.probe'.format(prefix, len(sampleID_LIST), len(passProbeID_DICT.keys()), call_rate_threshold * 100), 'w')
        for probeID in probeID_LIST:
            if probeID in passProbeID_DICT:
                fout.write(probeID + '\n')
        fout.close()

        fout = open('{0}_S{1:05}_P{2:05}_CR{3:02.0f}.geno'.format(prefix, len(sampleID_LIST), len(passProbeID_DICT.keys()), call_rate_threshold * 100), 'w')
        fout.write('\t'.join(['probeset_id'] + sampleID_LIST) + '\n')
        for probeID in probeID_LIST:
            if not probeID in passProbeID_DICT: continue
            context = [probeID]

            for sampleID in sampleID_LIST:
                context += [self.geno_DICT[sampleID][probeID]]

            fout.write('\t'.join(context) + '\n')
        fout.close()
        print('Write fileName: {0}, probeN: {1}, sampleN: {2}'.format(prefix, len(passProbeID_DICT.keys()), len(sampleID_LIST)))

    def make_hmp(self, prefix, probeInfo_DICT, probeID_LIST, call_rate_threshold=0, sampleID_LIST=None):
        if sampleID_LIST is None:
            sampleID_LIST = self.sampleID_LIST

        passProbeID_DICT = {}
        for probeID in probeID_LIST:
            if not probeID in self.probeID_DICT: continue
            if self.probeID_DICT[probeID] != self.fileN: continue

            noCallN = 0
            for sampleID in sampleID_LIST:
                if not probeID in self.geno_DICT[sampleID]:
                    noCallN += 1
                elif self.geno_DICT[sampleID][probeID] == 'NN':
                    noCallN += 1

            call_rate = 1 - (noCallN / len(sampleID_LIST))

            if call_rate >= call_rate_threshold:
                passProbeID_DICT[probeID] = call_rate

        fout = open('{0}_S{1:05}_P{2:05}_CR{3:02.0f}.probe'.format(prefix, len(sampleID_LIST), len(passProbeID_DICT.keys()), call_rate_threshold * 100), 'w')
        for probeID in probeID_LIST:
            if probeID in passProbeID_DICT:
                fout.write(probeID + '\n')
        fout.close()

        fout = open('{0}_S{1:05}_P{2:05}_CR{3:02.0f}.hmp.tab'.format(prefix, len(sampleID_LIST), len(passProbeID_DICT.keys()), call_rate_threshold * 100), 'w')
        fout.write('\t'.join(['rs', 'alleles', 'chrom', 'pos', 'strand', 'assembly', 'center', 'protLSID', 'assayLSID', 'panelLSID', 'QCcode'] + sampleID_LIST) + '\n')
        for probeID in probeID_LIST:
            if not probeID in passProbeID_DICT: continue

            if  probeInfo_DICT[probeID]['variantType'] == 'SNP' :
                alleles = probeInfo_DICT[probeID]['Allele_A'] + '/' + probeInfo_DICT[probeID]['Allele_B']
            else:
                alleles = 'A/T'
            chrom = probeInfo_DICT[probeID]['chrom']
            pos = probeInfo_DICT[probeID]['pos']
            strand = '+'
            #strand = probeInfo_DICT[probeID]['strand']
            assembly = 'NA'
            center = 'NA'
            protLSID = 'NA'
            assayLSID = 'NA'
            panelLSID = 'NA'
            QCcode = 'NA'

            context = [probeID, alleles, chrom, pos, strand, assembly, center, protLSID, assayLSID, panelLSID, QCcode]

            for sampleID in sampleID_LIST:
                context += [self.geno_DICT[sampleID][probeID]]

            fout.write('\t'.join(context) + '\n')
        fout.close()
        print('Write fileName: {0}, probeN: {1}, sampleN: {2}'.format(prefix, len(passProbeID_DICT.keys()), len(sampleID_LIST)))

#probeInfo_DICT
probeInfo_DICT = {}
fin = open('probe.info')
legend_LIST = fin.readline().rstrip('\n').split('\t')
for line in fin:
    data_LIST = line.rstrip('\n').split('\t')
    fields = {key: value for key, value in zip(legend_LIST, data_LIST)}
    probeID = fields['probeset_id']
    probeInfo_DICT[probeID] = fields
fin.close()


fin = open('probe.list')
probeID_LIST = []
probeID_DICT = {}
for line in fin:
    probeID = line.rstrip('\n')
    probeID_DICT[probeID] = 0
    probeID_LIST += [probeID]
fin.close()
print('probeN:', len(probeID_LIST))

fin = open('sample.list')
sampleID_LIST = []
sampleID_DICT = {}
for line in fin:
    sampleID = line.rstrip('\n')
    sampleID_DICT[sampleID] = 0
    sampleID_LIST += [sampleID]
fin.close()
print('sampleN:', len(sampleID_LIST))


geno = GENO()
geno.read_geno('0466')

#geno.make_geno('all', probeID_LIST, 0, sampleID_LIST)
#geno.make_fasta('all', probeID_LIST, 0, sampleID_LIST)
geno.make_hmp('all', probeInfo_DICT, probeID_LIST, 0, sampleID_LIST)
