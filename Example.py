from Pops import *


# Вид выражения, должен быть составлен по правилам синтаксиса

formula = [Sin()]
#formula = [Plus(), Sin(), PImp()]

# Параметры (для всех синусов и импульсов по 3 параметра по порядку нахождения в выражении)

params = [1, 1, 0]
#params = [0.1, 1, 0, 1, 10, 5]

# количество поколений
n_generation = 100

# количество особей в популяции
n_pops = 100

# количество мутаций
n_mut = 100

# количество скрещиваний
n_cross = 100

# количество лучших особей для показа
n_best = 3

# формирования ряда, который алгоритм будет пытаться угадать
y = Gen(formula)
y.put_gen_params(params)
y = y.value()


pops = Pops()
form = [Plus(), Sin(), PImp()]

for k in range(n_generation):
	pops.new_generation(n_pops=n_pops//2, chrom=formula)

	pops.crossover(n_cross)

	pops.mutate(n_mut)

	pops.fitsort(y)

	pops.leave_unique()

	pops.leave_strong(n_pops//2)


pops.leave_strong(n_best)

print('pars: ', params)
print('\n')
for i in pops.pops:
	i.obj.print_formula()
	print(i.fit)
	i.show(y)


