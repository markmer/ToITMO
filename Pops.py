from Gens import *





class Pop:

	def __init__(self, obj=Gen(), chrom=[] ,lifetime=0, fit=None):
		self.obj = obj
		self.chrom = chrom
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

	def fitness(self, y):
		if len(self.chrom) == 0:
			self.take_gen_params()
		self.put_gen_params()
		self.fit = np.linalg.norm(self.obj.value() - y)

	def mutation(self, eps=1):
		"""мутация с eps = 0 даст копию, мб сделать if с eps=None"""
		if self.fit == None:
			eps = eps
		else:
			eps = eps

		pop = copy.deepcopy(self.chrom)

		lp = len(pop)
		flag = True
		while flag:
			for i in range(lp):
				if r() == 0:
					pop[i] *= np.random.uniform(1-eps, 1+eps)
					flag = False

		new_pop = Pop(obj=self.obj.copy(), chrom=pop)
		new_pop.put_gen_params()
		return new_pop

	mutations = [mutation]


	def cross1(self, pop1):
		pop = copy.deepcopy(self.chrom)
		lp = len(pop)
		a = r(0, lp)
		pop[a] = 0.5*(pop[a] + pop1.chrom[a])

		new_pop = Pop(obj=self.obj.copy(), chrom=pop)
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

		new_pop = Pop(obj=self.obj.copy(), chrom=pop)
		new_pop.put_gen_params()
		return new_pop

	cross = [cross1, cross2]

	def show(self, y=0*t):
		self.obj.show(y)



class Pops:

	def __init__(self, pops=[]):
		self.pops = pops

	def new_random_generation(self, n_pops=1, n_gene=1):
		for i in range(n_pops):
			gen = Gen()
			gen.gene(n_gene)
			pop = Pop(gen)
			pop.init_params()
			self.pops += [pop]

	def new_generation(self, chrom=[Sin()], n_pops=1):
		gen = Gen(chrom)
		for i in range(n_pops):
			pop = Pop(gen.copy())
			pop.init_params()
			self.pops += [pop]

	def crossover(self, n_cross=1):
		pops1 = []
		lc = len(Pop.cross)
		lp = len(self.pops)
		for i in range(n_cross):
			r1 = r(0, lp)
			r2 = r(0, lp, r1)
			pops1 += [ Pop.cross[r(0, lc)](self.pops[r1], self.pops[r2]) ]
		self.pops = np.append(self.pops, pops1)

	def mutate(self, n_mut=1):
		pops1 = []
		lm = len(Pop.mutations)
		lp = len(self.pops)
		for i in range(n_mut):
			 pops1 += [ Pop.mutations[r(0, lm)](self.pops[r(0, lp)]) ]
		self.pops = np.append(self.pops, pops1)

	def sorting(self, l1):
		idx = np.argsort(l1)
		self.pops = [self.pops[i] for i in idx]

	def fitsort(self, y):
		fits = []
		for i in self.pops:
			if i.fit == None:
				i.fitness(y)
			fits += [i.fit]
		self.sorting(fits)

	def leave_strong(self, n_pops=1):
		self.pops = self.pops[:n_pops]

	def leave_unique(self):
		n = []
		uniq_pops = []
		flag = True
		for i in range(len(self.pops)):
			for j in n:
				if j == self.pops[i].fit:
					flag = False
					break
			if flag:
				n += [self.pops[i].fit]
				uniq_pops += [self.pops[i]]
			flag = True
		self.pops = uniq_pops





		
