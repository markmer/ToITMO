import numpy as np
from matplotlib.pyplot import *
import copy



t = np.arange(0, 100, 0.1)


def r(a=0, b=2, c=None):
	rn = np.random.randint(a, b)
	while rn == c:
		rn = np.random.randint(a, b)
	return rn

class Funcs:
	"""docstring for Funcs"""
	def __init__(self, p=[], const=False):
		self.p = np.array(p)
		self.const = const

	len_p = 1

	def func(self):
		return self.p

	def value(self):
		return self.func()
	
	def name(self):
		return type(self).__name__

	def init_params(self, a=0, b=1):
		self.p = np.random.uniform(a, b, self.len_p)

	def change_params(self, p):
		self.p = p[:self.len_p]

	def stay_const(self, flag=True):
		self.const = flag


class Sin(Funcs):
	"""docstring for Sin"""
	def __init__(self, p=[], const=False):
		super().__init__(p, const)

	len_p = 3

	def value(self):
		return self.p[0]*np.sin(self.p[1]*t + self.p[2])

	def init_params(self, aA=0, bA=10, aw=0, bw=0.2/abs(t[1]-t[0])):
		fi = np.random.uniform(0, 2*np.pi)
		w = np.random.uniform(aw, bw)
		A = np.random.uniform(aA, bA)
		self.p = np.array([A, w, fi])


class PImp(Funcs):
	"""docstring for PImp"""
	def __init__(self, p=[], const=False, symmetry=True):
		super().__init__(p, const)
		self.symmetry = symmetry

	len_p = 3

	def init_params(self, aA=0, bA=10, aT=t[1]-t[0], bT=t[-1]):
		T = np.random.uniform(aT, bT)
		tau = np.random.uniform(aT, T)
		A = np.random.uniform(aA, bA)
		self.p = np.array([A, T, tau])

	def value(self):
		m = np.zeros(len(t))
		cond = (t % self.p[1] < self.p[2])
		m[cond] = 1
		if self.symmetry:
			m[~cond] = -1
		return self.p[0]*m



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
		
		

class Gen:
	def __init__(self, chrom=[]):
		self.chrom = chrom

	f = Funcs.__subclasses__()
	g = GFuncs.__subclasses__()


	def gene(self, gen_len=1):
		# цепочка не может быть четной длины, четные значения gen_len -> gen_len+1
		# вероятность образования цепи уменьшается с длиной
		count = 0
		pop = []
		while count != -1:

			if len(pop)+count+1 < gen_len:
				if r() == 0:
					i = self.g[r(0, len(self.g))]()
					count += 1
				else:
					i = self.f[r(0, len(self.f))]()
					count -= 1
				
				pop += [i]
			
			else:
				pop += [self.f[r(0, len(self.f))]()]
				count -= 1
		self.chrom = pop


	def mark(s):
		count = 0 # индекс конца левого поддерева в строке
		i0 = len(s)
		if len(s) < 2:
			return [count, count]
		else:
			for i in range(len(s)):
				if issubclass(type(s[i]),Funcs):
					count -= 1
				else:
					count += 1
				if count == 0:
					i0 = min(i0, i)
				if count == -1:
					return [i0, i]
		

	def decprint1(s):
		dev = Gen.mark(s)[0]

		if issubclass(type(s[0]),Funcs):
			return s[0].name() + '(' + str(s[0].p) + ')'
		else:
			return s[0].name() + '(' + Gen.decprint1(s[1:dev+1]) + ',' + Gen.decprint1(s[dev+1:])  + ')'
	
	
	def print_formula(self):
		if len(self.chrom) == 0:
			print(None)
			return
		s = self.chrom[:]
		print(Gen.decprint1(s))
		

	def dec1(s):
		dev = Gen.mark(s)[0]

		if issubclass(type(s[0]),Funcs):
			return s[0].value()
		else:
			return s[0].value(Gen.dec1(s[1:dev+1]), Gen.dec1(s[dev+1:]))

	def value(self):
		if len(self.chrom) == 0:
			return None
		s = self.chrom[:]
		return Gen.dec1(s)

	def init_params(self):
		if len(self.chrom) == 0:
			return
		for i in self.chrom:
			i.init_params()

	def take_gen_params(self):
		par = []
		for i in self.chrom:
			if issubclass(type(i),Funcs) and not i.const:
				par += [i.p]
		return np.array(par).ravel()

	def put_gen_params(self, par):
		par = np.array(par)
		for i in self.chrom:
			if len(par) == 0:
				break
			if issubclass(type(i),Funcs) and not i.const:
				lp = len(par)
				if lp >= i.len_p:
					i.p = par[:i.len_p]
				else:
					i.p[:lp] = par[:lp]
				par = par[i.len_p:]

	def stay_const(self, flag=True):
		for i in self.chrom:
			i.stay_const(flag)

	def copy(self):
		return copy.deepcopy(self)

	def show(self, y=0*t):
		if len(self.chrom) == 0:
			print('empty Gen chrom for plot')
			return
		figure(str(self))
		plot(t, y)
		plot(t, self.value())
		show()

	