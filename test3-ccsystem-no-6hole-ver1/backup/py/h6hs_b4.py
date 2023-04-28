# import functools
import itertools

class SatProblem:
    def __init__(self):
        self.clauses = set()
        self.vars = {}
        self.n = 0
        self.known = {}
        self.comments = ['']

    def add_var(self, *ids):
        if ids not in self.vars:
            self.n += 1
            self.vars[ids] = self.n
        return self.vars[ids]
    
    def add_known(self, value, *ids):
        idx = self.add_var(*ids)
        self.set_known(value, idx)

    def set_known(self, value, idx):
        if idx < 0:
            idx, value = -idx, not value
        if idx not in self.known:
            self.known[idx] = value
        assert(self.known[idx] == value)
        
    def set_var(self, *ids):
        self.vars[ids[:-1]] = ids[-1]

    def var(self, *ids):
        return self.vars[ids]

    def add_clause(self, *ids):
        # self.clauses.add(tuple(sorted(ids)))
        nds = []
        value = False
        s = set(ids)
        for i in s:
            # if -i in s:
            #     print(i, ids)
            #     return 
            if i in self.known:
                value = value or self.known[i]
            elif -i in self.known:
                value = value or not self.known[-i]
            else:
                nds.append(i)
        if not value:
            assert(len(nds) > 0)
            self.clauses.add(tuple(sorted(nds)))

    def size(self):
        return self.n, len(self.clauses), sum(len(c) for c in self.clauses), len(self.known)

    def to_strings(self):
        for s in self.comments:
            yield 'c' + s
        n, m, _, k = self.size()
        yield 'p cnf %d %d' % (n - k, m)
        ids = {}
        def gget(id):
            if id == 0: return 0
            if abs(id) not in ids:
                idx = len(ids) + 1
                ids[abs(id)] = idx
            if id > 0:
                return ids[id]
            else:
                return -ids[-id]
        for c in self.clauses:
            yield ' '.join(map(lambda id: str(gget(id)), c)) + ' 0'

    def save(self, filename):
        with open(filename, 'w') as out:
            for s in self.to_strings():
                out.write(s + '\n')
            out.close()

    def define_var_and(self, idx, *vars):
        # self.add_clause(idx, *[-i for i in vars])
        # for i in vars:
        #     self.add_clause(-idx, i)
        value = True
        nds = []
        for i in vars:
            if i in self.known:
                value = value and self.known[i]
            elif -i in self.known:
                value = value and not self.known[-i]
            else:
                nds.append(i)
        if not value:
            self.set_known(False, idx)
        else:
            if len(nds) == 0:
                self.set_known(True, idx)
            else:
                self.add_clause(idx, *[-i for i in nds])
                for i in nds:
                    self.add_clause(-idx, i)

class HullStructure:
    def __init__(self, hs):
        self.n = sum(hs)
        self.m = len(hs)
        self.hs = hs
        self.color = []
        self.hss = []
        for i in range(len(hs)):
            self.hss.append([len(self.color) + j for j in range(hs[i])])
            self.color += [i] * hs[i]
        
class Summary:
    def __init__(self):
        self.s = ''
    def add(self, s):
        self.s += '\n' + s
    def __str__(self):
        return self.s

def ss_system(n, hull, smi):
    sat = SatProblem()
    for i, j, k in itertools.combinations(range(n), 3):
        idx = sat.add_var(i, j, k)
        sat.set_var(j, k, i, idx)
        sat.set_var(k, i, j, idx)
        sat.set_var(i, k, j, -idx)
        sat.set_var(k, j, i, -idx)
        sat.set_var(j, i, k, -idx)

    # convex
    for idx in range(hull.m):
        ids = hull.hss[idx]
        for i, j, k in itertools.combinations(ids, 3):
            sat.add_known(True, i, j, k)
        for i, j in zip(ids, ids[1:] + [ids[0]]):
            for k in range(ids[0], n):
                if k not in [i, j]:
                    sat.add_known(True, i, j, k)
    # symmetry
    for idx in range(hull.m - 1):
        ids = hull.hss[idx]
        nds = hull.hss[idx + 1]
        if len(ids) == 3:
            if len(nds) < 2:
                continue
            x, y = nds[:2]
            p, q = ids[:2]
            sat.add_known(True, x, p, y)
            sat.add_known(True, x, y, q)
    smi.c3v, _, _, smi.c3vk = sat.size()

    for p, q, r, t in itertools.permutations(range(n), 4):
        sat.add_clause(-sat.var(t, q, r), -sat.var(p, t, r), -sat.var(p, q, t), sat.var(p, q, r))
    for p, q, r, s, t in itertools.permutations(range(n), 5):
        sat.add_clause(-sat.var(t, s, p), -sat.var(t, s, q), -sat.var(t, s, r), -sat.var(t, p, q), -sat.var(t, q, r), sat.var(t, p, r))
    return sat

def define_pt_inside_triangle(sat, n, hull):
    for p, q, r, s in itertools.permutations(range(n), 4):
        if p > q or p > r:
            continue
        idx = sat.add_var(p, q, r, s)
        sat.set_var(q, r, p, s, idx)
        sat.set_var(r, p, q, s, idx)
        colors = [hull.color[i] for i in [p, q, r, s]]
        cp, cq, cr, cs = colors
        if cs <= min(cp, cq, cr):
            sat.set_known(False, idx)
        else:
            sat.define_var_and(idx, sat.var(p, q, r), sat.var(p, q, s), sat.var(q, r, s), sat.var(r, p, s))

def add_clasuse_hole(sat, n, pts):
    vars = []
    for i in range(len(pts)):
        p, q = pts[i - 1], pts[i]
        for r in pts:
            if r != p and r != q:
                vars.append(-sat.var(p, q, r))
    p, q, r, s = pts
    #   s--r  
    # L |  | R
    #   p--q  
    lnopts, rnopts = [], []
    for x in range(n):
        if x in pts:
            continue
        # for i in range(2, len(pts) - 2):
        for i in range(1, len(pts) - 1):
            vars.append(sat.var(pts[0], pts[i], pts[i + 1], x))
    
        lidx = sat.add_var(p, q, r, s, x, -1) # left
        sat.define_var_and(lidx, sat.var(p, q, x), sat.var(r, s, x), sat.var(p, s, x))
        lnopts.append(-lidx)
        
        ridx = sat.add_var(p, q, r, s, x, -2) # right
        sat.define_var_and(ridx, sat.var(p, q, x), sat.var(r, s, x), sat.var(r, q, x))
        rnopts.append(-ridx)

    lpts_idx = sat.add_var(p, q, r, s, -1)
    sat.define_var_and(lpts_idx, *lnopts)
    rpts_idx = sat.add_var(p, q, r, s, -2)
    sat.define_var_and(rpts_idx, *rnopts)
    vars += [lpts_idx, rpts_idx]
    sat.add_clause(*vars)

def add_clause_all_hole(sat, n, m):
    for p in range(n):
        for pts in itertools.permutations(range(p + 1, n), m - 1):
            add_clasuse_hole(sat, n, [p] + list(pts))

def all_hull_structure(n):
    def hs_dfs(n, now):
        if n == 0 and now[-1] < 6:
            yield now
        if n < 3:
            if n == 1 and now[-1] < 8 or n == 2 and now[-1] < 7:
                yield now + [n]
        for m in range(3, min(n + 1, 9)):
            for hs in hs_dfs(n - m, now + [m]):
                yield hs
    for hs in hs_dfs(n, []):
        yield hs



n = 14

print("n =", n)
hss = set()
for hs in all_hull_structure(n):
    hss.add(tuple(hs))
    if len(hss) == 1:
        print(hs)
        summary = Summary()
        h = HullStructure(hs)
        sat = ss_system(n, h, summary)
        sat.comments.append('hs: ' + ' '.join(map(str, hs)))
        print(n, 'size ss', sat.size())
        define_pt_inside_triangle(sat, n, h)
        print(n, 'size tr', sat.size())
        add_clause_all_hole(sat, n, 4)
        print(n, 'size tt', sat.size())
        filename = 'no_6hole_n%d_h%d_b4.sat' % (n, len(hss))
        sat.save(filename)
print(len(hss))
