import sys


nb_dimensions=[5, 10,15];
fonctions = ['Ackley ', 'Quadric ', 'Griewank ', 'Sphere ', 'Rastrigin ','Rosenbrock ', 'Schwefel1_2 ','Schwefel ', 'Salomon ' ];

f = open('task_gwo.jobs','w')

for j in nb_dimensions:
	for fonction in fonctions:
		f.write('./exemplecluster '+ fonction +str(j)+' spso2011 \n' )

f.write('exit\n')
f.close()
