from gen import *


class GenAlg:
	def __init__(self, krit_fit=0.5, krit_lifetime=100):
		self.gens = []
		self.bank = []
		self.krit_fit = krit_fit
		self.krit_lifetime = krit_lifetime

	def add_gen(self, gen):
		self.gens += [gen]

	def new_random_generation(self, n_gens=1, n_gene=1):
		for i in range(n_gens):
			gen = Gen()
			gen.gene(n_gene)
			gen.init_params()
			self.gens += [gen]

	def increase_lifetime(self):
		for i in self.gens:
			i.increase_lifetime()

	def sorting(self, l1, for_bank=False):
		idx = np.argsort(l1)
		if for_bank:
			self.bank = [self.bank[i] for i in idx]
		else:
			self.gens = [self.gens[i] for i in idx]

	def fitsort(self, y, norm_y=1, norm_len=7, norm_worth= 0.9):
		fits = []
		for i in self.gens:
			if i.fit == None:
				i.fitness(y=y, norm_y=norm_y, norm_len=norm_len, norm_worth=norm_worth)
			fits += [i.fit]
		self.sorting(fits)

	def bank_fitsort(self):
		fits = []
		for i in self.bank:
			fits += [i.fit]
		self.sorting(l1=fits, for_bank=True)

	def crossover(self, n_cross=1, best=False):
		new_gens = []
		lc = len(self.gens[0].crosss)
		lp = len(self.gens)
		if best:
			for i in range(n_cross):
				r1 = r(1, lp)
				new_gens += [ self.gens[0].crosss[r(0, lc)](self.gens[r1], self.gens[0]) ]
		else:
			for i in range(n_cross):
				r1 = r(0, lp)
				r2 = r(0, lp, r1)
				new_gens += [ self.gens[0].crosss[r(0, lc)](self.gens[r1], self.gens[r2]) ]
		self.gens = np.append(self.gens, new_gens)

	def mutate(self, n_mut=1, best=False):
		new_gens = []
		lm = len(self.gens[0].mutations)
		lp = len(self.gens)
		if best:
			for i in range(n_mut):
				 new_gens += [ self.gens[0].mutations[r(0, lm)](self.gens[0]) ]
		else:
			for i in range(n_mut):
				 new_gens += [ self.gens[0].mutations[r(0, lm)](self.gens[r(0, lp)]) ]
		self.gens = np.append(self.gens, new_gens)

	def kill_best_area_fit(self, krit=0.01):
		best_gen = self.gens[0]
		for i in self.gens[1:]:
			if abs(i.fit - best_gen.fit) <= krit:
				if i.lifetime > best_gen.lifetime:
					best_gen.lifetime = i.lifetime
				self.gens.remove(i)

	def kill_long_gens(self, maxx=7):
		a = []
		for i in self.gens:
			if len(i.chrom) <= maxx:
				a += [i]
		self.gens = a

	def kill_all_with_weak_best_fit(self, a=2):
		best_gen = self.gens[0]

		flag1 = (best_gen.fit > self.krit_fit) and (best_gen.lifetime > self.krit_lifetime)
		flag2 = best_gen.lifetime > a*self.krit_lifetime

		if flag1 or flag2:
			if flag1:
				print()
				print('flag1')
				best_gen.print_formula()
				print(best_gen.fit, best_gen.lifetime)
				print('Популяция была убита')

				if len(self.bank) == 0:
					self.krit_fit = best_gen.fit
				self.bank += [best_gen]
				self.leave_strong(0)

				print('new krit_lifetime', self.krit_lifetime)
				print('new krit_fit', self.krit_fit)
				print()
				return

			relative_fit = self.krit_fit/best_gen.fit

			if best_gen.lifetime > a*relative_fit*self.krit_lifetime:
				print()
				print('----> flag2 <----')
				best_gen.print_formula()
				print(best_gen.fit, best_gen.lifetime)
				print('Популяция была убита')

				self.bank += [best_gen]
				self.krit_lifetime *= relative_fit
				self.krit_fit = best_gen.fit
				self.leave_strong(0)

				print('new krit_lifetime', self.krit_lifetime)
				print('new krit_fit', self.krit_fit)
				print()
				return

	def leave_strong(self, n_gens=1):
		self.gens = self.gens[:n_gens]

	def leave_unique(self):
		n = []
		uniq_gens = []
		flag = True
		for i in range(len(self.gens)):
			for j in n:
				if j == self.gens[i].fit:
					flag = False
					break
			if flag:
				n += [self.gens[i].fit]
				uniq_gens += [self.gens[i]]
			flag = True
		self.gens = uniq_gens

	def stop_GA(self, stop=5, eps_stop=0.05):
		if len(self.bank) >= stop:
			self.bank_fitsort()
			count = 0
			for i in self.bank:
				if abs(i.fit - self.bank[0].fit) <= eps_stop:
					count += 1
				else:
					break
			if count >= stop:
				return True
			else:
				return False
		else:
			return False

	
	"""----------------------------------------------------------"""
	def GA(self, y, norm_y=1, n_generation=100, n_gens=100, n_cross=20, n_mut=20, n_strong = 50, maxlen=10):

		self.leave_strong(0)
		for k in range(n_generation):
			self.new_random_generation(n_gens=n_gens - len(self.gens), n_gene=3)
			self.crossover(n_cross)
			self.mutate(n_mut)
			self.kill_long_gens(maxlen)
			self.fitsort(y=y, norm_y=norm_y, norm_len=maxlen, norm_worth=0.9)
			self.leave_unique()
			self.leave_strong(n_strong)
			self.increase_lifetime()

			flag = False
			if self.gens[0].fit < 0.02:
				flag = True
			if flag:
				break

		print('количество итераций', k)
		for i in self.gens[:2]:
			i.print_formula()
			print(i.fit, i.lifetime)

		if self.gens[0].fit < 0.1:
			print('yes', self.gens[0].fit)
		else:
			print('no', self.gens[0].fit)
		print('\n')

		return self.gens[0]


	def GA_best(self, y, norm_y=1, n_generation=100, n_gens=100, n_cross=20, n_mut=20, n_strong = 50, \
		best_area_fit=0.01, stop=5, stop_eps=0.005, stop_fit=0.05, maxlen=7):
		
		self.new_random_generation(n_gens=n_gens - len(self.gens))
		self.fitsort(y=y, norm_y=norm_y)
		flag = False
		for k in range(n_generation):
			self.crossover(n_cross, True)
			self.mutate(n_mut, True)
			self.kill_long_gens(maxlen)
			self.fitsort(y=y, norm_y=norm_y, norm_len=maxlen, norm_worth=0.9)
			self.kill_best_area_fit(best_area_fit)
			self.leave_strong(n_strong)
			self.increase_lifetime()

			if len(self.bank) != 0 and len(self.bank)%20 == 0:
				self.krit_lifetime *= 1.2

			flag = self.stop_GA(stop, stop_eps)

			if self.gens[0].fit < stop_fit or flag:
				if not flag:
					self.bank += [self.gens[0]]
				break

			self.kill_all_with_weak_best_fit()
			self.new_random_generation(n_gens=n_gens - len(self.gens))


		self.bank += [self.gens[0]]
		print('количество итераций', k)
		print('В банке : ', len(self.bank), 'особей')
		self.bank_fitsort()
		for i in self.bank:
			i.print_formula()
			print(i.fit, i.lifetime)
		print('\n')

		return self.bank[0]




class GenAlgPop(GenAlg):
	"""-----------------------"""
	def __init__(self, pop_pattern=Gen(), krit_fit=0.5, krit_lifetime=100):
		super().__init__(krit_fit, krit_lifetime)
		self.pop_pattern = pop_pattern

	def new_random_generation(self, n_gens=1, n_gene=1):
		for i in range(n_gens):
			gen = self.pop_pattern.copy()
			gen.init_params()
			self.gens += [gen]

	def kill_best_area_norm(self, krit=1):
		best_pop = self.pops[0]
		for i in self.pops[1:]:
			if np.linalg.norm(i.chrom - best_pop.chrom) < krit:
				self.pops.remove(i)

	def kill_bank_area_norm(self, krit=1):
		if len(self.bank) == 0:
			return
		else:
			for i in self.pops:
				for j in self.bank:
					if np.linalg.norm(i.chrom - j.chrom) < krit:
						self.pops.remove(i)
	def kill_long_gens(self):
		pass

	def GA_best(self, y, norm_y=1, n_generation=100, n_gens=15, n_cross=3, n_mut=3, n_strong = 7, best_area_fit=0.01,\
	 stop=5, stop_eps=0.001, stop_fit=0.05):
		
		n_params = self.pop_pattern.obj.count_params()
		n_gens *= n_params
		n_cross *= n_params
		n_mut *= n_params
		n_strong *= n_params

		self.new_random_generation(n_gens=n_gens - len(self.gens))
		self.fitsort(y=y, norm_y=norm_y)
		flag = False
		for k in range(n_generation):
			self.crossover(n_cross, True)
			self.mutate(n_mut, True)
			self.kill_long_gens()
			self.fitsort(y=y, norm_y=norm_y)
			self.kill_best_area_fit(best_area_fit)
			self.leave_strong(n_strong)
			self.increase_lifetime()

			if len(self.bank) != 0 and len(self.bank)%20 == 0:
				self.krit_lifetime *= 1.2

			flag = self.stop_GA(stop, stop_eps)

			if self.gens[0].fit < stop_fit or flag:
				if not flag:
					self.bank += [self.gens[0]]
				break

			self.kill_all_with_weak_best_fit()
			self.new_random_generation(n_gens=n_gens - len(self.gens))


		self.bank += [self.gens[0]]
		print('количество итераций', k)
		print('В банке : ', len(self.bank), 'особей')
		self.bank_fitsort()
		for i in self.bank:
			i.obj.print_formula()
			print(i.fit, i.lifetime)
		print('\n')

		return self.bank[0]