from optparse import OptionParser
import sys
import gzip
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'input',help = "")
parser.add_option("-p","--prefix",action = 'store',type = 'string',dest = 'prefix',help = "")
parser.add_option("-n","--readN",action = 'store',type = 'int',dest = 'readN',help = "")

(opt, args) = parser.parse_args()
if opt.input == None or opt.prefix == None or opt.readN == None:
    print('Basic usage')
    print('')
    print('     python split_fastq.py -i test_1.fastq -p test_1.fastq -n 100000')
    print('')
    sys.exit()

#infile = 'Keumgang.01cell.OmniC_1.fastq.gz'
#readN = 100000
infile = opt.input
readN  = opt.readN
prefix = opt.prefix

fin = gzip.open(infile, 'rt')

fileN = 0

fout = None
for lineIDX, line in enumerate(fin):
    if lineIDX%(readN*4) == 0:
        fileN += 1
        if fout != None: fout.close()
        fout = gzip.open(prefix + '.{:08}'.format(fileN), 'wt')
    fout.write(line)
fout.close()
