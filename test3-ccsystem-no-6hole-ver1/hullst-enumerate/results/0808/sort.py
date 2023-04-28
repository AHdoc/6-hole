f = open("unsat_s_0808.txt", "r")

s = []
for line in f.readlines():
    for x in line.split():
        s.append(x)

s.sort()
print(s)
