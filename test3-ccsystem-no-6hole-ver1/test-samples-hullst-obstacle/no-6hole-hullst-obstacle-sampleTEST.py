import subprocess
import time

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def test(nb, n, s):
    start = time.time()
    result_generate = subprocess.getoutput("./a-hullst-obstacle "+s).split()
    result_sat = subprocess.getoutput("time ./cadical-sc2020-45029f8/build/cadical -n -q no-6hole-"+s+"-obstacle.sat").split()
    result_remove = subprocess.getoutput("rm no-6hole-"+s+"-obstacle.sat").split()
    end = time.time()

    out = str(nb)+" "+str(n)+" "+s+": "+result_sat[1]+" "+str(end - start)+' '+time.ctime(end)
    if result_sat[1] == "UNSATISFIABLE":
        print(bcolors.WARNING + out + bcolors.ENDC)
        return False
    else:
        print(out)
        return True

f = open("n_sample.txt", "r")
for line_f in f.readlines():
    tmp = line_f.split()
    nb = int(tmp[0])
    n = int(tmp[1])
    s = tmp[2]
    test(nb, n, s)
