f_uns_old = open("unsat_s_old.txt","r")
g = open("n_le_24.txt","r")
h_sat = open("sat_s.txt","w")
h_uns = open("unsat_s.txt","w")

uns_old = []
for line in f_uns_old.readlines():
    for s in line.split():
        uns_old.append(s)

sat = []
uns = []

for line in g.readlines():
    tmp = line.split()
    s = tmp[2]

    chk_old = True
    for s2 in uns_old:
        if s2 in s:
            chk_old = False
    
    if chk_old:
        sat.append(s)
    if (not chk_old) and (s in uns_old) and s.endswith('0'):
        uns.append(s)

for s in sat:
    h_sat.write(s+' ')
h_sat.write('\n')

for s in uns:
    h_uns.write(s+' ')
h_uns.write('\n')

f_uns_old.close()
g.close()
h_sat.close()
h_uns.close()

