def check(f_filename, g_filename):
    f = open(f_filename, "r")
    g = open(g_filename, "r")

    unsat_s = []
    for line in f.readlines():
        for s in line.split():
            unsat_s.append(s)

    for line in g.readlines():
        for s in line.split():
            chk = True
            for s2 in unsat_s:
                if s2 in s:
                    chk = False
            if not chk:
                print('ERROR '+s)
            #else:
            #    print('OK '+s)

print('...')
check("unsat_s.txt", "25.txt")
print('---')
