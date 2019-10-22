import numpy as np
from matplotlib.pyplot import *
import copy



t = np.arange(0, 100, 0.1)


def r(a=0, b=2, c=None):
	rn = np.random.randint(a, b)
	while rn == c:
		rn = np.random.randint(a, b)
	return rn

def geom_r(c, p=0.05):
	rn = np.random.geometric(p)
	while rn >= c:
		rn = np.random.geometric(p)
	return rn




class Funcs:
	"""docstring for Funcs"""
	def __init__(self, p=[], const=False):
		self.p = np.array(p)
		self.const = const

	def period(self, f, t, T):
		for i in range(int(t/T)):
			if i == int(t/T) - 1:
				f = np.append(f, f[:t % T])
			else:
				f = np.append(f, f[:T])
		return f
		
	def name(self):
		return type(self).__name__

	def change_params(self, p):
		self.p = p[:self.len_p]

	def stay_const(self, flag=True):
		self.const = flag

	def copy(self):
		return copy.deepcopy(self)


class Sin(Funcs):
	"""docstring for Sin"""
	def __init__(self, p=[], const=False):
		super().__init__(p, const)

		self.len_p = 3

	def value(self):
		self.p = np.array(self.p, dtype=float)
		if self.p[0] < 0:
			self.p[0] = -self.p[0]
			self.p[2] = self.p[2] + np.pi

		if self.p[1] < 0:
			self.p[1] = -self.p[1]
			self.p[2] = -self.p[2] + np.pi

		if abs(self.p[2]) > 2*np.pi:
			self.p[2] = self.p[2] - abs(2*np.pi)*np.sign(self.p[2])*(abs(self.p[2]) // abs(2*np.pi))
		if self.p[2] < 0:
			self.p[2] = self.p[2] + 2*np.pi

		return self.p[0]*np.sin(self.p[1]*t + self.p[2])

	def init_params(self, aA=0, bA=10, aw=0, bw=0.2/abs(t[1]-t[0])):
		if not self.const: # and len(self.p) < self.len_p:
			fi = np.random.uniform(0, 2*np.pi)
			w = np.random.uniform(aw, bw)
			A = np.random.uniform(aA, bA)
			self.p = np.array([A, w, fi])


class PImp(Funcs):
	"""docstring for PImp"""
	def __init__(self, p=[], const=False, symmetry=True):
		super().__init__(p, const)
		self.symmetry = symmetry

		self.len_p = 3

	def init_params(self, start_A=0, stop_A=10, start_T=t[1]-t[0], stop_T=t[-1]):
		if not self.const:
			T = np.random.uniform(start_T, stop_T)
			tau = np.random.uniform(start_T, T)
			A = np.random.uniform(start_A, stop_A)
			self.p = np.array([A, T, tau])

	def value(self):
		if self.p[2] > self.p[1]:
			self.p[2] = self.p[1]

		m = np.zeros(len(t))
		cond = (t % self.p[1] < self.p[2])
		m[cond] = 1
		if self.symmetry:
			m[~cond] = -1
		return self.p[0]*m


class TImp(Funcs):
	"""--------------"""
	def __init__(self, p=[], const=False):
		super().__init__(p, const)

		self.len_p = 4

	def init_params(self, start_A=0, stop_A=10, start_T=t[1]-t[0], stop_T=t[-1]):
		if not self.const:
			T = np.random.uniform(start_T, stop_T)
			tau = np.random.uniform(start_T, T)
			t0 = np.random.uniform(start_T, tau)
			A = np.random.uniform(start_A, stop_A)
			self.p = np.array([A, T, tau, t0])

	def value(self):
		if self.p[2] > self.p[1]:
			self.p[2] = self.p[1]
		if self.p[3] > self.p[2]:
			self.p[3] = self.p[2]

		t1 = t[t <= self.p[1]]
		cond1 = (t1 < self.p[3])
		cond2 = (t1 >= self.p[3]) & (t1 <= self.p[2])
		m = np.zeros(len(t1))
		m[cond1] = t1[cond1]/self.p[3]
		if self.p[3] == self.p[2]:
			m[cond2] = 1 + self.p[3]/(self.p[2] + 1-self.p[3]) - t1[cond2]/(self.p[2] + 1 -self.p[3])
		else:
			m[cond2] = 1 + self.p[3]/(self.p[2]-self.p[3]) - t1[cond2]/(self.p[2]-self.p[3])
		return self.p[0]*self.period(m, len(t), len(t1))






class GFuncs:
	"""docstring for GFuncs"""
	def __init__(self):
		pass

	def name(self):
		return type(self).__name__

	def init_params(self):
		pass

	def stay_const(self, flag=True):
		pass

	def copy(self):
		return copy.deepcopy(self)

class Plus(GFuncs):
	"""docstring for Plus"""
	def __init__(self):
		pass

	def value(self, *f):
		return f[0]+f[1]

class Product(GFuncs):
	def __init__(self):
		pass

	def value(self, *f):
		return f[0]*f[1]
