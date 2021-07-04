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

unsat_s = []
sat_s = []

def test(nb, n, s):
    start = time.time()
    result_generate = subprocess.getoutput("./a-hullst-obstacle "+s).split()
    result_sat = subprocess.getoutput("time ./cadical-sc2020-45029f8/build/cadical -n -q no-6hole-"+s+"-obstacle.sat").split()
    result_remove = subprocess.getoutput("rm no-6hole-"+s+"-obstacle.sat").split()
    end = time.time()

    out = str(nb)+" "+str(n)+" "+s+": "+result_sat[1]+" "+str(end - start)+' '+time.ctime(end)
    if result_sat[1] == "UNSATISFIABLE":
        print(bcolors.WARNING + out + bcolors.ENDC)
        global unsat_s
        unsat_s.append(s)
        return False
    else:
        print(out)
        global sat_s
        sat_s.append(s)
        return True

def solve(nb, n, s):
    global unsat_s
    global sat_s

    unsat_s.clear()
    g = open("unsat_s.txt", "r")
    for line_g in g.readlines():
        for s2 in line_g.split():
            unsat_s.append(s2)
    g.close()

    sat_s.clear()
    g = open("sat_s.txt", "r")
    for line_g in g.readlines():
        for s2 in line_g.split():
            sat_s.append(s2)
    g.close()

    chk = False
    for s2 in unsat_s:
        if s2 in s:
            chk = True
            break
    for s2 in sat_s:
        if s in s2:
            chk = True
            break

    if not chk:
        if test(nb, n, s):
            h = open("sat_s.txt", "w")
            for s2 in sat_s:
                h.write(s2+' ')
            h.write('\n')
            h.close()
        else:
            h = open("unsat_s.txt", "w")
            for s2 in unsat_s:
                h.write(s2+' ')
            h.write('\n')
            h.close()


#idx = 0
#for s in ['877710', '777630', '777540', '7774410', '677730', '577740', '5377710', '477750', '4777410', '4477710', '3777510', '3577710', '33844710']:
#    idx += 1
#    solve(idx, 30, s)

f = open("n_21_to_24.txt", "r")
for line_f in f.readlines():
    tmp = line_f.split()
    nb = int(tmp[0])
    n = int(tmp[1])
    s = tmp[2]

    solve(nb, n, s)




