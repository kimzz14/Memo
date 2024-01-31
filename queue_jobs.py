from optparse import OptionParser
import subprocess,time,os,shlex,smtplib,sys,time
from datetime import timedelta
devnull = open(os.devnull,'w')

def refresh(proc_LIST):
    _process = subprocess.Popen('ps -o pid,command', stdout=subprocess.PIPE, shell=True)
    out, err = _process.communicate()
    line_LIST = out.decode('utf-8').strip('\n').split('\n')
    status_DICT = {}
    for line in line_LIST[1:]:
        line = line.lstrip(' ')
        pid = line.split(' ')[0]
        command = line.lstrip(pid + ' ')
        status_DICT[int(pid)] = command
    for proc in proc_LIST:
        if proc.status == 'running':
            if proc.pid not in status_DICT: 
                proc.status = 'complete'
            elif status_DICT[proc.pid].find('defunct') != -1: 
                proc.process.kill()
                proc.status = 'complete'

class PROC:
    def __init__(self, command, pid, status):
        self.process = None
        self.command = command
        self.pid = pid
        self.status = status

    def run(self):
        global devnull
        self.process = subprocess.Popen(self.command, stdout=devnull, shell=True)
        self.pid = self.process.pid
        self.status = 'running'

    def refresh(self):
        if self.status == 'complete':
            pass
        elif self.status == 'ready':
            pass
        elif self.status == 'running':
            _process = subprocess.Popen('ps --pid ' +  str(self.pid) + ' -o command', stdout=subprocess.PIPE, shell=True)
            out, err = _process.communicate()
            command = out.decode('utf-8').split('\n')[1]
            if command == '':
                self.status = 'complete'
            elif command.find('defunct') != -1:
                #self.process.kill()
                self.status = 'complete'
        else:
            print('bug!')
        return self.status

parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-p","--pid",action = 'store',type = 'string',dest = 'PID',help = "")
parser.add_option("-m","--message",action = 'store',type = 'string',dest = 'MESSAGE',help = "")
parser.add_option("-i","--input",action = 'store',type = 'string',dest = 'INPUT',help = "")
parser.add_option("-n","--maxProcN",action = 'store',type = 'int',dest = 'maxProcN',help = "")
parser.add_option("-t","--prefix",action = 'store',type = 'string',dest = 'PREFIX',help = "")
(opt, args) = parser.parse_args()
if opt.INPUT == None or opt.maxProcN == None or opt.PREFIX == None:
    print('Basic usage')
    print('')
    print("     queue_jobs.py --input command.sh --prefix NGS01 --maxProcN 10 --pid '1953,1235' ")
    print('')
    sys.exit()

if opt.PID is None:
    prePID_LIST = []
else:
    prePID_LIST = list(map(int, opt.PID.strip(' ').rstrip('\n').rstrip(',').split(',')))
inFile = opt.INPUT
prefix = opt.PREFIX
maxProcN = opt.maxProcN


readyIDX = 0
allProcN_LIST = []
runProc_LIST = []
readyProc_LIST = []
completeProc_LIST = []

#prePID
for prePID in prePID_LIST:
    proc = PROC('NA', prePID, 'running')
    runProc_LIST += [proc]

#add Command
fin = open(inFile)
for line in fin:
    if line.startswith('#') == True: continue
    command = line.rstrip('\n')
    proc = PROC(command, None, 'ready')
    readyProc_LIST += [proc]
    allProcN_LIST  += [proc]
fin.close()


#fout_log
fout_log = open(prefix + '.log','w')

tic = time.time()
while True:
    isChange = False
    #refresh
    refresh(runProc_LIST)
    #print('a', len(runProc_LIST))
    _runProc_LIST = []
    for procIDX, proc in enumerate(runProc_LIST):
        if proc.status == 'running': 
            _runProc_LIST += [proc]
        elif proc.status == 'complete': 
            completeProc_LIST += [proc]
            isChange = True

    runProc_LIST = _runProc_LIST
    #print('b', len(runProc_LIST))

    toc = time.time()
    processTime = toc - tic
    remainProcN = len(allProcN_LIST) - len(completeProc_LIST)

    if len(completeProc_LIST) == 0: 
        remainProcessTime = 0
    else:
        remainProcessTime = (processTime / len(completeProc_LIST)) * remainProcN


    isDone = False
    for idx in range(maxProcN - len(runProc_LIST)):
        proc = readyProc_LIST[readyIDX]
        runProc_LIST += [proc]
        proc.run()
        readyIDX += 1
        if len(readyProc_LIST) == readyIDX:
            isDone = True
            break

    #print('c', len(runProc_LIST))
    if isChange:
        #save command
        context  = []
        context += [time.strftime("%H:%M:%S")]
        context += [str(timedelta(seconds=int(processTime))).rjust(20, ' ')]
        context += [str(timedelta(seconds=int(remainProcessTime))).rjust(20, ' ')]
        context += [(str(len(completeProc_LIST)) + '/' + str(len(allProcN_LIST))).rjust(14, ' ')]
        context = '  '.join(map(str,context))
        fout_log.write(context + '\n')
        fout_log.flush()
    if isDone == True:
        break
    time.sleep(1)
fout_log.close()
