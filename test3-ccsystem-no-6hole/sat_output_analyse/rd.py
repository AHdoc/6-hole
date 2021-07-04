f = open('no-6hole-4466710-obstacle-1points-extra-simplified.sat')

line = f.readline()
while line:
    if line.startswith('379 ') or line.startswith('-379 ') or ' 379 ' in line or ' -379 ' in line:
        print(line)
        #break
    line = f.readline()
