from fractions import Fraction
from functools import cmp_to_key
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.collections import PatchCollection
import time
import random
import cProfile
import math


def sign(x):
	if x == 0: return 0
	return 1 if x > 0 else -1

class Point:
	def __init__(self, x, y):
		self.x = Fraction(x)
		self.y = Fraction(y)
	def __neg__(self):
		return Point(-self.x, -self.y)
	def __add__(self, other):
		return Point(self.x + other.x, self.y + other.y)
	def __sub__(self, other):
		return Point(self.x - other.x, self.y - other.y)
	def __mul__(self, other):
		return Point(self.x * other, self.y * other)
	def __truediv__(self, other):
		return Point(self.x / other, self.y / other)
	def __hash__(self):
		return hash((self.x, self.y))
	def __str__(self):
		return "[%s, %s]" % (str(self.x), str(self.y))
	def __lt__(self, other):
		if self.x == other.x:
			return self.y < other.y
		return self.x < other.x 
	def __eq__(self, other):
		return self.x == other.x and self.y == other.y
	def __le__(self, other):
		return self < other or self == other
	def Dot(self, other):
		return self.x * other.x + self.y * other.y
	def Det(self, other):
		return self.x * other.y - self.y * other.x
	def Rot90(self):
		return Point(-self.y, self.x)
	def Len2(self):
		return Point.Dot(self, self)
	def Length(self):
		return self.Len2() ** 0.5
	def Unit(self):
		return self / Fraction(self.Length())

	@staticmethod
	def Cross(a, b, c):
		return Point.Det(b - a, c - a)
	@staticmethod
	def CrossOp(a, b, c):
		return sign(Point.Cross(a, b, c))
	@staticmethod
	def CounterClockwise(a, b):
		A = a.x > 0 or a.x == 0 and a.y > 0
		B = b.x > 0 or b.x == 0 and b.y > 0
		if A != B:
			return 1 if A else -1
		d = a.Det(b)
		if d == 0:
			if a == b:
				return 0
			return 1 if a.Len2() > b.Len2() else -1
		return 1 if d > 0 else -1

	@staticmethod
	def StrList(a):
		return "[%s]"% (", ".join([str(p) for p in a]))
	@staticmethod
	def Rectangle(x1, y1, x2, y2):
		return [Point(x1, y1), Point(x2, y1), Point(x2, y2), Point(x1, y2)]
	@staticmethod
	def GetRange(ps):
		xs, ys = [p.x for p in ps], [p.y for p in ps]
		return [min(xs), min(ys), max(xs), max(ys)]
	
	@staticmethod
	def PlotList(a, sty):
		plt.plot([p.x for p in a], [p.y for p in a], sty)
	@staticmethod
	def PlotPolygons(polygons):
		fig, ax = plt.subplots()
		patches = [matplotlib.patches.Polygon([[p.x, p.y] for p in polygon]) for polygon in polygons]
		# print(len(patches))
		pc = PatchCollection(patches, alpha=1.0, facecolors=("black",))
		# pc.set_array(np.array([0] * len(polygons)))
		ax.add_collection(pc)

def ConvexHull(ps): # 求点集 ps 组成的凸包
	n = len(ps)
	if n <= 1: return ps
	ps = sorted(ps)
	# print(Point.StrList(ps))
	qs = []
	for p in ps:
		while len(qs) > 1 and Point.CrossOp(qs[-2], qs[-1], p) <= 0: qs.pop()
		qs.append(p)
	t = len(qs)
	for p in ps[-2:0:-1] + [ps[0]]:
		while len(qs) > t and Point.CrossOp(qs[-2], qs[-1], p) <= 0: qs.pop()
		qs.append(p)
	qs.pop()
	return qs

def GetIntersectionPoint(p1, p2, q1, q2):
	a1, a2 = Point.Cross(q1, q2, p1), -Point.Cross(q1, q2, p2)
	# print(Point.StrList([p1, p2, q1, q2]), a1, a2)
	return (p1 * a2 + p2 * a1) / (a1 + a2)

def IsParallel(p1, p2, q1, q2):
	return sign((p2 - p1).Det(q2 - q1)) == 0

def IsColinear(u, v, p):
	return Point.CrossOp(p, u, v) == 0

def IsOnSegment(u, v, p): 
	return Point.CrossOp(p, u, v) == 0 and sign(Point.Dot(u - p, v - p)) <= 0


def ConvexCut(ps, a, b): #多边形半平面交
	qs = []
	for i in range(len(ps)):
		p, q = ps[i], ps[(i+1)%len(ps)]
		if Point.Cross(a, b, p) > 0 or IsColinear(a, b, p):
			qs.append(p)
		if not IsParallel(p, q, a, b):
			t = GetIntersectionPoint(p, q, a, b)
			if IsOnSegment(p, q, t):
				qs.append(t)
	return ConvexHull(qs)

def IsContain(polygon, p):
	for i in range(len(polygon)):
		x = polygon[i]
		y = polygon[i + 1] if i + 1 < len(polygon) else polygon[0]
		if Point.Cross(x, y, p) < 0:
			return False
	return True

def CheckInclusion(polygon, ps):
	for p in ps:
		if p not in polygon and IsContain(polygon, p):
			return True
	return False

def CheckValidComplete(ps):
	n = len(ps)
	print(n, "points to check")
	for i1 in range(n):
		for i2 in range(i1):
			if ps[i1] == ps[i2]:
				print("Same points! @ ", i1, i2, ps[i1], ps[i2])
				return False
	print("All points distinct")
	for i1 in range(n):
		for i2 in range(i1):
			for i3 in range(i2):
				if Point.CrossOp(ps[i1], ps[i2], ps[i3]) == 0:
					print("Colinear points! @ ", i1, i2, i3, ps[i1], ps[i2], ps[i3])
					return False
	print("No colinear points")
	for i1 in range(n):
		for i2 in range(i1):
			for i3 in range(i2):
				for i4 in range(i3):
					convex = ConvexHull([ps[i1], ps[i2], ps[i3], ps[i4]])
					if len(convex) < 4 or CheckInclusion(convex, ps):
						continue
					for i5 in range(i4):
						convex = ConvexHull([ps[i1], ps[i2], ps[i3], ps[i4], ps[i5]])
						if len(convex) < 5 or CheckInclusion(convex, ps):
							continue
						for i6 in range(i5):
							ids = [i1, i2, i3, i4, i5, i6]
							convex = ConvexHull([ps[i] for i in ids])
							if len(convex) == 6 and not CheckInclusion(convex, ps):
								print("Found empty convex hexagon @ ", ids, Point.StrList(convex))
								return False
	print("No empty convex hexagon exists")
	return True

def GetInit(n = 29):
	raw = [ [ 1,1260], [ 16, 743], [ 22, 531], [ 37, 0], [ 306, 592], [ 310, 531], [ 366, 552], [ 371, 487], [ 374, 525], [ 392, 575], [ 396, 613], [ 410, 539], [ 416, 550], [ 426, 526], [ 434, 552], [ 436, 535], [ 446, 565], [ 449, 518], [ 450, 498], [ 453, 542], [ 458, 526], [ 489, 537], [ 492, 502], [ 496, 579], [ 516, 467], [ 552, 502], [ 754, 697], [ 777, 194], [1259, 320] ]
	# raw = raw[:-1]
	# raw = [ [0, 0], [1, 1] ]
	ret = [Point(a[0], a[1]) for a in raw]
	return ret

def GetGeneralPositions(n, limit = 100, seed = 0):
	random.seed(seed)
	ps = []
	for i in range(n):
		while True:
			p = Point(random.randrange(0, limit), random.randrange(-limit, limit))
			ok, invalid = True, False
			for i1 in range(len(ps)):
				ok = ok and (p != ps[i1] and p.x != ps[i1].x)
				for i2 in range(i1):
					ok = ok and not IsColinear(ps[i1], ps[i2], p)
					for i3 in range(len(ps)):
						ok = ok and not IsParallel(ps[i1], ps[i2], ps[i3], p)
						if i3 < i2 and len(ConvexHull([p, ps[i1], ps[i2], ps[i3]])) == 4:
							for i4 in range(i3):
								for i5 in range(i4):
									convex = ConvexHull([p, ps[i1], ps[i2], ps[i3], ps[i4], ps[i5]])
									if len(convex) == 6 and not CheckInclusion(convex, ps):
										invalid = True
			if ok and not invalid:
				ps.append(p)
				break
		print(len(ps), ps[-1])
	return ps

def GetGravityCenter(polygon):
	ss, sa = Point(0, 0), 0
	for i in range(len(polygon)):
		x, y = polygon[i], polygon[i + 1] if i + 1 < len(polygon) else polygon[0]
		a = x.Det(y)
		ss += (x + y) / 3 * a
		sa += a
	return ss / sa

def GetSimpleCenter(polygon):
	rs = [p.x.denominator * p.y.denominator for p in polygon]
	ss, sa = Point(0, 0), 0
	for i in range(len(polygon)):
		ss += polygon[i] * rs[i]
		sa += rs[i]
	return ss / sa

def SamplePointsInPolygon(polygon, k = 1):
	ret = [GetSimpleCenter(polygon)]
	while len(ret) < k:
		tmp = []
		for g in ret:
			for i in range(len(polygon)):
				x, y = polygon[i], polygon[i + 1] if i + 1 < len(polygon) else polygon[0]
				tmp.append((g + x + y) / 3)
		ret += tmp
	if len(ret) > k:
		ret = ret[:k]
	return ret

def GetLastH5s(ps, p):
	ret = []
	for i1 in range(len(ps)):
		for i2 in range(i1):
			for i3 in range(i2):
				C4 = ConvexHull([p, ps[i1], ps[i2], ps[i3]])
				if len(C4) == 4:
					for i4 in range(i3):
						C5 = ConvexHull([p, ps[i1], ps[i2], ps[i3], ps[i4]])
						if len(C5) == 5 and not CheckInclusion(C5, ps):
							ret.append(C5)
	return ret

def GetAllH5s(ps):
	ret = []
	for i in range(len(ps)):
		ret += GetLastH5s(ps[:i], ps[i])
	return ret

limit = 10 ** 8
limit_range = Point.Rectangle(-limit, -limit, limit, limit)

def GetFarPoint(A, B, R = limit_range):
	ret, d = A, 0
	for i in range(len(R)):
		x, y = R[i], R[i + 1] if i + 1 < len(R) else R[0]
		if not IsParallel(A, B, x, y):
			o = GetIntersectionPoint(A, B, x, y)
			od = (o - A).Dot(B - A)
			if od > 0 and (d == 0 or od < d):
				ret, d = o, od
	return ret

def GetSingleBlackArea(A, B, C, D):
	#     E.F
	#    /   \
	#   C-----B
	#  /       \
	# D         A
	if not IsParallel(A, B, C, D):
		o = GetIntersectionPoint(A, B, C, D)
		if (o - A).Dot(B - A) > 0:
			R = limit_range
			if IsContain(R, o):
				return [B, o, C]
	# E = C - (D - C).Unit() * 100
	# F = B + (B - A).Unit() * 100
	E, F = GetFarPoint(D, C), GetFarPoint(A, B)
	return [B, F, E, C]

def GetBlackAreas(h5s):
	# assert(ps == sorted(ps))
	blacks = []
	for ids in h5s:
		tmp = ids + ids
		for i in range(len(ids)):
			A, B, C, D = tmp[i : i + 4]
			blacks.append(GetSingleBlackArea(A, B, C, D))
	return blacks


def GetBlackAreasStrict(h5s):
	# assert(ps == sorted(ps))
	blacks = []
	def get(A, B, C, D):
		#     E.F
		#    /   \
		#   C-----B
		#  /       \
		# D         A
		if not IsParallel(A, B, C, D):
			o = GetIntersectionPoint(A, B, C, D)
			if (o - A).Dot(B - A) > 0:
				R = limit_range
				if IsContain(R, o):
					return [B, o, C]
		# E = C - (D - C).Unit() * 100
		# F = B + (B - A).Unit() * 100
		E, F = GetFarPoint(D, C), GetFarPoint(A, B)
		return [B, F, E, C]

	for ids in h5s:
		tmp = ids + ids
		for i in range(len(ids)):
			A, B, C, D = tmp[i : i + 4]
			blacks.append(GetSingleBlackArea(A, B, C, D))
	return blacks

def GetHalfLine(A, B, C, D):
	# get the left part of line A, B on line C, D
	P = GetIntersectionPoint(A, B, C, D)
	s = 1
	if (B - A).Det(D - C) < 0:
		s = -1
	return P, s

def MergeHalfLines(a, D):
	L, R = None, None
	for P, s in a:
		if s > 0:
			if L is None or (P - L).Dot(D) > 0:
				L = P
		else:
			if R is None or (P - R).Dot(D) < 0:
				R = P
	if L is not None and R is not None and (L - R).Dot(D) > 0:
		return None, None
	return L, R

class PointPool:
	def __init__(self):
		self.list = []
		self.dict = {}
	def Get(self, x):
		if x not in self.dict:
			self.dict[x] = (len(self.list), [])
			self.list.append(x)
		return self.dict[x]
	def Index(self, x):
		return self.Get(x)[0]
	def Set(self, x, y):
		itr = self.Get(x)
		self.dict[x] = (itr[0], y)
	def __getitem__(self, x):
		return self.Get(x)[1]

class BlackArea:
	def __init__(self, A, B, C, D):
		self.A = A
		self.B = B
		self.C = C
		self.D = D
		self.points = [A, B, C, D]
		self.o = GetIntersectionPoint(A, B, C, D)
		self.open = (self.o - A).Dot(B - A) < 0
	def Right(self, x):
		return self.open or self.o.x > x
	def GetPolygon(self, x = None, limit = limit):
		Range = Point.Rectangle(-limit, -limit, limit, limit)
		L, R = self.B, self.C
		X, Y = Point(x, 0), Point(x, 1)
		if x is not None:
			L = GetIntersectionPoint(self.A, self.B, X, Y)
			R = GetIntersectionPoint(self.C, self.D, X, Y)
		if not self.open and IsContain(Range, self.o):
			return [L, self.o, R]
		# E = C - (D - C).Unit() * 100
		# F = B + (B - A).Unit() * 100
		E, F = GetFarPoint(self.D, self.C, Range), GetFarPoint(self.A, self.B, Range)
		if x is not None and self.C < self.D:
			R = GetFarPoint(X, Y, Range)
			return [R, F, L]
		if x is not None and self.B < self.A:
			L = GetFarPoint(Y, X, Range)
			return [L, E, R]
		return [L, F, E, R]
	def IsContain(self, p):
		return Point.Cross(self.A, self.B, p) >= 0 and Point.Cross(self.C, self.B, p) >= 0 and Point.Cross(self.C, self.D, p) >= 0
		
def GetBlackRight(convex, x):
	'''Get restrains on the right side of rightmost point and x'''
	i = (convex.index(max(convex)) - 2) % len(convex)
	A, B, C, D, E = convex[i:] + convex[:i]
	ret = [BlackArea(A, B, C, D), BlackArea(B, C, D, E)]
	ret = [b for b in ret if b.Right(x)]
	return ret

def DebugPlots(ps, points, blacks, infpoints, edges = []):
	polygons = [b.GetPolygon(x = ps[-1].x) for b in blacks]
	Range = Point.GetRange(points.list)
	polygons.append(Point.Rectangle(Range[0], Range[1], ps[-1].x, Range[3]))
	polygons.append(Point.Rectangle(Range[0], Range[1], ps[-1].x, Range[3]))
	Point.PlotPolygons(polygons)
	Point.PlotList(ps, 'ro')
	Point.PlotList([p for p in points.list if p not in infpoints], 'rx')
	Point.PlotList([p for p in points.list if p in infpoints], 'go')
	for X, Y in edges:
		D = (Y - X).Rot90().Unit() / 10
		Point.PlotList([X + D, Y + D], '')
	plt.xlim(-10, 2000)
	plt.ylim(-2000, 2000)

class Line:
	def __init__(self, i1, i2, P, Q):
		self.i1 = i1
		self.i2 = i2
		self.P = P
		self.Q = Q
		self.qs = []
		self.qsid = []
		self.D = Q - P

def SearchNewPointRaw(ps):
	if len(ps) <= 2:
		raise Exception('point set size %d too small' % len(ps))
	points = PointPool()
	Z, ZY = ps[-1], ps[-1] + Point(0, 1)
	for i1 in range(len(ps)):
		points.Index(ps[i1])
	lines = []
	for i1 in range(len(ps)):
		for i2 in range(i1):
			lines.append(Line(i1, i2, ps[i2], ps[i1]))
	lines.append(Line(len(ps) - 1, len(ps) - 1, Z, ZY))
	infpoints = set()
	for line in lines:
		for line_other in lines:
			if not IsParallel(line.Q, line.P, line_other.P, line_other.Q):
				line.qs.append(GetIntersectionPoint(line.Q, line.P, line_other.P, line_other.Q))
		line.qs.sort()
		line.qs += [line.qs[-1] + line.qs[-1] - line.qs[0], line.qs[0] - line.qs[-1] + line.qs[0]]
		infpoints.add(line.qs[-1])
		infpoints.add(line.qs[-2])
		line.qs.sort()
	blacks = []
	for h5s in GetAllH5s(ps):
		blacks += GetBlackRight(h5s, Z.x)
	print("Pool size", len(points.list), 'Black count', len(blacks))
	# DebugPlots(ps, points, blacks, infpoints)
	# plt.show()

	edges = set()
	for line in lines:
		E = line.D.Rot90() / 10
		F = GetIntersectionPoint(line.P, line.Q, Z, ZY) if not IsParallel(line.P, line.Q, Z, ZY) else Z
		covers = []
		for B in blacks:
			L, R, left, right = None, None, True, True
			constrains = [(B.A, B.B), (B.C, B.B), (B.C, B.D), (Z, Z - Point(0, 1))]
			parallels = [IsParallel(x, y, line.P, line.Q) for x, y in constrains]
			halflines = [GetHalfLine(constrains[i][0], constrains[i][1], line.P, line.Q) for i in range(len(constrains)) if not parallels[i]]
			L, R = MergeHalfLines(halflines, line.Q - line.P)	
			if L is not None or R is not None:
				if True in parallels:
					A, B = constrains[parallels.index(True)]
					if (line.Q - line.P).Dot(B - A) > 0:
						right = False
					else:
						left = False
				if R is None: R = line.qs[-1]
				if L is None: L = line.qs[0]
				if R.x >= Z.x:
					if L.x < Z.x: L = F
					if L < R: covers.append((L, R, left, right))
		if line.P.x == line.Q.x:
			covers.append((line.qs[0], line.qs[-1], True, False))
		covers.sort()
		nqsid, nqs = [], []
		cid = 0
		ur, dr = None, None
		for i in range(len(line.qs)):
			if Z.x <= line.qs[i].x and (len(nqs) == 0 or nqs[-1] != line.qs[i]):
				nqsid.append(points.Index(line.qs[i]))
				nqs.append(line.qs[i])
				if len(nqsid) > 1:
					if ur is None:
						edges.add((nqsid[-2], nqsid[-1]))
					if dr is None:
						edges.add((nqsid[-1], nqsid[-2]))
				if ur is not None and ur <= line.qs[i]: ur = None
				if dr is not None and dr <= line.qs[i]: dr = None
				while cid < len(covers) and covers[cid][0] <= line.qs[i]:
					if covers[cid][2] and (ur is None or ur < covers[cid][1]):
						ur = covers[cid][1]
					if covers[cid][3] and (dr is None or dr < covers[cid][1]):
						dr = covers[cid][1]
					cid += 1
		line.qs, line.qsid = nqs, nqsid	
	g = [[] for i in range(len(points.list))]
	
	# DebugPlots(ps, points, blacks, infpoints)
	# for x, y in edges:
	# 	X, Y = points.list[x], points.list[y]
	# 	D = (Y - X).Rot90().Unit() / 10
	# 	Point.PlotList([X + D, Y + D], '')
	# plt.show()
	
	next_edge = {}
	for i, j in edges:
		g[i].append((j, False))
		g[j].append((i, True))
	for i in range(len(g)):
		if len(g[i]) > 0:
			g[i].sort(key = lambda p: (cmp_to_key(Point.CounterClockwise)(points.list[p[0]] - points.list[i]), p[1]))
			g[i].append(g[i][0])
			for j in range(1, len(g[i])):
				x, y = g[i][j - 1], g[i][j]
				if x[1] and not y[1]:
					next_edge[(x[0], i)] = y[0]
	# extract blocks
	blocks = []
	visit = set()
	for a, b in next_edge:
		if points.list[a] in infpoints:
			ids = [a, b]
			na, nb = a, b
			visit.add((a, b))
			while points.list[nb] not in infpoints:
				na, nb = nb, next_edge[(na, nb)]
				visit.add((na, nb))
				ids.append(nb)
			blocks.append(ids)
	infblkcnt = len(blocks)
	
	for a, b in next_edge:
		if (a, b) not in visit:
			ids = [a]
			na, nb = a, b
			while (na, nb) not in visit:
				visit.add((na, nb))
				ids.append(nb)
				na, nb = nb, next_edge[(na, nb)]
			ids = ids[:-1]
			blocks.append(ids)
	ret = []
	for block in blocks:
		polygon = [points.list[i] for i in block]
		ret += SamplePointsInPolygon(polygon, 1)

	polygons = [b.GetPolygon(x = Z.x) for b in blacks]
	polygons = polygons + polygons
	Range = Point.GetRange(points.list)
	polygons.append(Point.Rectangle(Range[0], Range[1], Z.x, Range[3]))
	polygons.append(Point.Rectangle(Range[0], Range[1], Z.x, Range[3]))
	for block in blocks:
		tmp = [points.list[i] for i in block]
		G = GetGravityCenter(tmp)
		polygons.append([p + (G - p).Unit() for p in tmp])
	Point.PlotPolygons(polygons)
	Point.PlotList(ps, 'ro')
	Point.PlotList(points.list, 'rx')
	Point.PlotList([p for p in points.list if p in infpoints], 'go')
	plt.xlim(-10, 5000)
	plt.ylim(-5000, 5000)
	plt.show()
	
	# print("%d blocks in total, of which %d blocks infinite" % (len(blocks), infblkcnt))
	return ret


stack = []
best = [[], 0]
def SearchNewPoint(depth):
	ps, _ = stack[depth]
	if len(ps) > len(best[0]):
		best[0], best[1] = ps, 1
		print('new size', len(ps), ':', Point.StrList(ps))
	# elif len(ps) == len(best[0]):
	# 	best[1] += 1
		# print('multiple size', len(ps), best[1], ':', Point.StrList(ps))
	candidates = SearchNewPointRaw(ps)
	for p in candidates:
		if len(stack) < depth + 2:
			stack.append([])
		stack[depth + 1] = [ps + [p], depth]
		SearchNewPoint(depth + 1)		


def TestBlack(s = ''):
	# ps = GetInit()
	ps = []
	for t in s.replace(')', ' ').replace('(', ' ').strip().split():
		x, y = map(int, t.split(','))
		ps.append(Point(x, y))
	# ps = GetGeneralPositions(6)
	# ps = [Point(a[0], a[1]) for a in [[0, 0], [3, 0], [4, 2], [1, 5], [0, 4]]]
	# print(Point.StrList(Point.Rectangle(-limit, -limit, limit, limit)))
	ps = sorted(ps)
	h5s = GetAllH5s(ps)
	print(len(h5s), "5-holes")
	blacks = GetBlackAreas(h5s)
	print(len(blacks), "black areas")
	Point.PlotPolygons(blacks)
	Point.PlotList(ps, 'ro')
	plt.show()

# def TestSimple():
# 	ps = GetInit()[:17]
# 	# ps = GetGeneralPositions(14, 1000)
# 	# print(CheckValidComplete(ps))
# 	ps = [Point(a[0], a[1]) for a in [[-1, 4], [0, 0], [1, 5], [2, 12], [3, 0], [4, 3]]]
# 	ps = sorted(ps)
# 	SearchNewPointRaw(ps)


# cProfile.run('TestSimple()', sort='cumtime')

ff = open('log_first_29.txt').readlines()
# ff = open('tmp.txt').readlines()
sss = []
for si in range(len(ff)):
	ss = ff[si].strip()
	if ss.startswith('29:'):# and int(ss[:-1])>=23:
		s = ff[si][3:].strip()
		sss.append((si, s))
#import random
#random.shuffle(sss)
for si in range(len(sss)):
	s = sss[si][1]
	osi = sss[si][0]
	# TestSimple()
	# s = '(67272,86876)(-12196,123572)(18522,10892)(0,0)(-53292,-53846)(-98150,-188154)(-14529,-7632)(-591,8014)(59879,88862)(5003,-13204)(33570,1033)(-12107,-13472)(-3937,18884)(114758,171976)(191720,35526)(-61802,-33373)(-217544,46081)(-47331,-125880)(152521,159089)(198197,204547)(8265,-10710)(-72637,-154216)(-124414,-237227)(-45016,14522)(-30893,222132)'
	ps = []
	for t in s.replace(')', ' ').replace('(', ' ').strip().split():
		x, y = map(int, t.split(','))
		ps.append(Point(x, y))
	ps = sorted(ps)
	ps0 = ps[:]
	print(osi, si, len(ps), Point.StrList(ps))
	# ps = ps[:-1]
	# ps = []
	# for i in range(5):
	# 	ps.append(Point( int(math.cos(i*math.pi*2/5)*100), int(math.sin(i*math.pi*2/5)*100)  ))
	# ps.append(Point(-48, 128))
	# ps = GetGeneralPositions(6)
	# ps = [Point(a[0], a[1]) for a in [[0, 0], [3, 0], [4, 2], [1, 5], [0, 4]]]
	# print(Point.StrList(Point.Rectangle(-limit, -limit, limit, limit)))
	c4s = []
	for i1 in range(len(ps)):
		for i2 in range(i1):
			for i3 in range(i2):
				for i4 in range(i3):
					C4 = ConvexHull([ps[i1], ps[i2], ps[i3], ps[i4]])
					if len(C4) == 4:
						c4s.append(C4)
	print('c4s', len(c4s))
	b4s = []
	for pp in c4s:
		# a, b, c, d = pp
		n = len(pp)
		empty = True
		ok = [False] * n
		# a --- d
		# |     |
		# |     |
		# b --- c
		for p in ps:
			if p in pp:
				continue
			cs = [(Point.Cross(pp[i], pp[(i + 1)%n], p) > 0) for i in range(n)]
			# print(Point.StrList(pp), cs)
			if all(cs):
				empty = False
				break
			for i in range(n):
				if not cs[i] and cs[(i+n-1)%n] and cs[(i+1)%n]:
					ok[i] = True
		# print('!!', empty, ok)
		if not empty:
			continue
		for i in range(n):
			if ok[i]:
				b4s.append((pp+pp)[i+1:i+1+n])

	print('b4s', len(b4s))

	def drawbs(A, B, C, D):
		for p in ps:
			if p in [A, B, C, D]:
				continue
			if Point.Cross(A, B, p) > 0 and Point.Cross(B, C, p) < 0 and Point.Cross(C, D, p) > 0:
				return drawbs(A, B, p, C) + drawbs(B, p, C, D)

		if not IsParallel(A, B, C, D):
			o = GetIntersectionPoint(A, B, C, D)
			if (o - A).Dot(B - A) > 0:
				R = limit_range
				if IsContain(R, o):
					return [[B, o, C]]
		# E = C - (D - C).Unit() * 100
		# F = B + (B - A).Unit() * 100
		# E, F = GetFarPoint(D, C), GetFarPoint(A, B)
		R = limit_range
		t = ConvexCut(R, A, B)
		t = ConvexCut(t, C, B)
		t = ConvexCut(t, C, D)
		return [t]
	blacks = []
	for bb in b4s:
		#     E.F
		#    /   \
		#   C-----B
		#  /       \
		# D         A
		a, b, c, d = bb
		blacks += drawbs(a, b, c, d)
					# for i4 in range(i3):
					# 	C5 = ConvexHull([p, ps[i1], ps[i2], ps[i3], ps[i4]])
					# 	if len(C5) == 5 and not CheckInclusion(C5, ps):
					# 		ret.append(C5)

	print('blacks', len(blacks))
	# h5s = GetAllH5s(ps)
	# print(len(h5s), "5-holes")
	# blacks = GetBlackAreas(h5s)
	# print(len(blacks), "black areas")
	Point.PlotPolygons(blacks)
	Point.PlotList(ps, 'ro')
	# Point.PlotList(ps0[-1:], 'bx')
	plt.show()

# TestBlack('(67272,86876)(-12196,123572)(18522,10892)(0,0)(-53292,-53846)(-98150,-188154)(-14529,-7632)(-591,8014)(59879,88862)(5003,-13204)(33570,1033)(-12107,-13472)(-3937,18884)(114758,171976)(191720,35526)(-61802,-33373)(-217544,46081)(-47331,-125880)(152521,159089)(198197,204547)(8265,-10710)(-72637,-154216)(-124414,-237227)(-45016,14522)(-30893,222132)')
# now = GetGeneralPositions(15, 100000)
# now = GetInit()
# now = sorted(now)
# now = now[:28]
# stack.append([now, 0])
# SearchNewPoint(0)
# now = sorted(GetInit())[:27]
# stack = []
# stack.append([now, 0])
# SearchNewPoint(0)
