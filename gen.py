from funcs import *


class Gen:
	def __init__(self, chrom=[], params=[], lifetime=0, fit=None):
		self.chrom = chrom
		self.params = params
		self.lifetime = lifetime
		self.fit = fit

		self.f = Funcs.__subclasses__()
		self.g = GFuncs.__subclasses__()
	
	def mark(s):
		count = 0 
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

	def take_gen_params(self):
		self.params = []
		for i in self.chrom:
			if issubclass(type(i),Funcs) and not i.const:
				for j in i.p:
					self.params += [j]
		self.params = np.array(self.params, dtype=float)
		return copy.deepcopy(self.params)

	def put_gen_params(self, par=[]):
		if len(par) == 0:
			pars = np.array(copy.deepcopy(self.params), dtype=float)
		else:
			pars = par
		for i in self.chrom:
			lp = len(pars)
			if lp == 0:
				break
			if issubclass(type(i),Funcs) and not i.const:
				if lp >= i.len_p:
					i.p = pars[:i.len_p]
				else:
					i.p[:lp] = pars[:lp]
				pars = pars[i.len_p:]

	def init_params(self):
		if len(self.chrom) == 0:
			return
		for i in self.chrom:
			i.init_params()
		self.take_gen_params()

	def count_params(self):
		count = 0
		for i in self.chrom:
			if issubclass(type(i),Funcs) and not i.const:
				count += i.len_p
		return count

	def increase_lifetime(self, n=1):
		self.lifetime += n

	def stay_const_all(self, flag=True):
		for i in self.chrom:
			i.stay_const(flag)

	def copy(self):
		return copy.deepcopy(self)

	def formprint(s):
		dev = Gen.mark(s)[0]

		if issubclass(type(s[0]),Funcs):
			return s[0].name() + '(' + str(s[0].p) + ')'
		else:
			return s[0].name() + '(' + Gen.formprint(s[1:dev+1]) + ',' + Gen.formprint(s[dev+1:])  + ')'
	
	
	def print_formula(self):
		if len(self.chrom) == 0:
			print(None)
			return
		s = self.chrom[:]
		print(Gen.formprint(s))

	def val(s):
		dev = Gen.mark(s)[0]

		if issubclass(type(s[0]),Funcs):
			return s[0].value()
		else:
			return s[0].value(Gen.val(s[1:dev+1]), Gen.val(s[dev+1:]))

	def value(self):
		if len(self.chrom) == 0:
			return None
		s = self.chrom[:]
		return Gen.val(s)

	def fitness(self, y, norm_y=1, norm_len=7, norm_worth= 0.9):
		if len(self.params) == 0:
			self.take_gen_params()
		self.put_gen_params()
		self.fit = norm_worth*np.linalg.norm(self.value() - y)/norm_y + (1 - norm_worth)*len(self.chrom)/norm_len

	def show(self, y=0*t):
		if len(self.chrom) == 0:
			print('empty Gen chrom for plot')
			return
		figure(str(self))
		plot(t, y)
		plot(t, self.value())
		show()


	# GA for chrom
	def chrom_names(self):
		names = []
		for i in self.chrom:
			names += [i.name()]
		return np.array(names)

	def cross(self, gen1):
		if len(self.chrom) == 1 and len(gen1.chrom) == 1:
			return self.append(gen1=gen1)

		if len(self.chrom) == 1:
			start = 0
		else:
			start = r(1, len(self.chrom))

		subtree = self.chrom[start:]
		stop = Gen.mark(subtree)[1]
		subtree = self.chrom[start:start+stop+1]

		start1 = r(0, len(gen1.chrom))
		subtree1 = gen1.chrom[start1:]
		stop1 = Gen.mark(subtree1)[1]
		subtree1 = gen1.chrom[start1:start1+stop1+1]

		new_chrom = np.append(np.array(self.chrom[:start]), np.array(subtree1))
		new_chrom = np.append(new_chrom, np.array(self.chrom[start+stop+1:]))

		new_gen = Gen(new_chrom)
		new_gen.take_gen_params()
		return new_gen.copy()		

	def append(self, gen1, g=None):
		new_chrom = np.append(np.array(self.chrom), np.array(gen1.chrom))
		if g == None:
			g = self.g[r(0, len(self.g))]()
		new_chrom = np.append(np.array(g), new_chrom)
		
		new_gen = Gen(new_chrom)
		new_gen.take_gen_params()
		return new_gen.copy()

	crosss = [cross, cross, cross, append]

	def mutation(self):
		new_gen_chrom = copy.deepcopy(self.chrom)
		# while Gen.formprint(new_gen_chrom) == Gen.formprint(self.chrom):
		for i in range(len(new_gen_chrom)):
			if r() == 0:
				if issubclass(type(new_gen_chrom[i]),Funcs):
					j = self.f[r(0, len(self.f))]()
					# while j.name() == new_gen_chrom[i].name():
					# 	j = self.f[r(0, len(self.f))]()
					j.init_params()
					new_gen_chrom[i] = j
				else:
					j = self.g[r(0, len(self.g))]()
					while j.name() == new_gen_chrom[i].name():
						j = self.g[r(0, len(self.g))]()
					new_gen_chrom[i] = j
		new_gen = Gen(new_gen_chrom)
		new_gen.take_gen_params()
		return new_gen

	# GA for params не создает копию, влияет на self
	def mutation_params(self, eps=1):
		if self.fit == None:
			eps = eps
		else:
			eps = self.fit
		
		lp = len(self.params)
		self.params = np.array(self.params, dtype=float)
		flag = True
		while flag:
			for i in range(lp):
				if r() == 0:
					self.params[i] *= np.random.uniform(1-eps, 1+eps)
					flag = False

	def mutation_full(self,eps=1):
		new_gen = self.mutation()
		new_gen.mutation_params(eps)
		new_gen.put_gen_params()
		return new_gen



	mutations = [mutation_full]



class Pop:

	def __init__(self, obj=Gen(), lifetime=0, fit=None):
		self.obj = obj
		self.chrom = []
		self.lifetime = lifetime
		self.fit = fit

	def take_gen_params(self):
		self.chrom = self.obj.take_gen_params()

	def put_gen_params(self):
		self.obj.put_gen_params(self.chrom)

	def init_params(self):
		self.obj.init_params()
		self.take_gen_params()

	def increase_lifetime(self, n=1):
		self.lifetime += n

	def put_in_chrom(self, chrom):
		self.chrom = chrom

	def print_formula(self):
		self.obj.print_formula()

	def fitness(self, y, norm_y=1, norm_len=7, norm_worth=0.9):
		if len(self.chrom) == 0:
			self.take_gen_params()
		self.put_gen_params()
		self.fit = np.linalg.norm(self.obj.value() - y)/norm_y

	def copy(self):
		return copy.deepcopy(self)

	def mutation(self, eps=1, fit0=10):
		if self.fit == None:
			eps = eps
		else:
			eps = self.fit

		pop = copy.deepcopy(self.chrom)

		lp = len(pop)
		flag = True
		while flag:
			for i in range(lp):
				if r() == 0:
					pop[i] *= np.random.uniform(1-eps, 1+eps)
					flag = False

		new_pop = Pop(obj=self.obj.copy())
		new_pop.chrom = pop
		new_pop.put_gen_params()
		return new_pop

	mutations = [mutation]

	def cross1(self, pop1):
		pop = copy.deepcopy(self.chrom)
		lp = len(pop)
		flag1 = True
		while flag1:
			for i in range(lp):
				if r() == 0:
					pop[i] = 0.5*(pop[i] + pop1.chrom[i])
					flag1 = False

		new_pop = Pop(obj=self.obj.copy())
		new_pop.chrom = pop
		new_pop.put_gen_params()
		return new_pop

	def cross2(self, pop1):
		pop = copy.deepcopy(self.chrom)
		lp = len(pop)
		flag1 = True
		flag2 = True
		while flag1 or flag2:
			for i in range(lp):
				if r() == 0:
					pop[i] = self.chrom[i]
					flag1 = False
				else:
					pop[i] = pop1.chrom[i]
					flag2 = False

		new_pop = Pop(obj=self.obj.copy())
		new_pop.chrom = pop
		new_pop.put_gen_params()
		return new_pop

	crosss = [cross1, cross2]

	def show(self, y=0*t):
		self.obj.show(y)