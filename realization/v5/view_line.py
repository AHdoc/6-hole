from fractions import Fraction
from functools import cmp_to_key
import matplotlib.pyplot as plt
import matplotlib.patches
from matplotlib.collections import PatchCollection
import time
import random
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

import os
if not os.path.exists('out'):
    os.makedirs('out')
ff = open('tmp_r8_shake_02081824.txt').readlines()
rs, ns = [], []
for si in range(len(ff)):
	ss = ff[si].strip()
	if ss.endswith(':'):
		s = ff[si + 1].strip()
		ps = []
		for t in s.replace(')', ' ').replace('(', ' ').strip().split():
			x, y = map(int, t.split(','))
			ps.append(Point(x, y))
		rs.append(len(rs))
		ns.append(len(ps))
		# ps = sorted(ps)
		# Point.PlotList(ps, 'ro')
		# fn = 'out/no%03d_%d.png'%(si, len(ps))
		# plt.xlim([-10000*30, 10000*30])
		# plt.ylim([-10000*30, 10000*30])
		# plt.savefig(fn)
		# print(si, len(ps), fn)
		# plt.clf()
plt.plot(rs, ns)
plt.show()
