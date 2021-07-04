def check(f_filename, g_filename):
    f = open(f_filename, "r")
    g = open(g_filename, "r")

    unsat_s = []
    for line in f.readlines():
        for s in line.split():
            unsat_s.append(s)

    cnt = [0 for i in range(10)]
    for line in g.readlines():
        tmp = line.split()
        s = tmp[2]
    
        chk = True
        for s2 in unsat_s:
            if s2 in s:
                chk = False
        
        if chk:
            cnt[ord(s[0])-ord('0')] += 1
            #if '25' in g_filename and s.startswith('6'):
            #    print(s)
    
    for i in range(3,9):
        print(str(i)+':'+str(cnt[i]), end='   ')
    print('   sum='+str(sum(cnt)))

for n in range(1,31):
    print('n='+str(n), end='   ')
    check("unsat_s.txt", "n_"+str(n)+"_smallest_hulls.txt")
