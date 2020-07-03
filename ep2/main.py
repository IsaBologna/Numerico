import numpy as np
import sys
import math

from ep_methods.tarefas import Teste, Teste_B
import ep_methods.numerical_methods as nm

N=128
T=1
p=np.array([0.35])

# batata = Teste(N,T,p)

# nm.matrix_uk(batata)
# nm.solve_normal_system(batata)
# print(batata.a)

p2=np.array([0.15,0.3,0.7,0.8])

batata2 = Teste_B(N,T,p2)

nm.matrix_uk(batata2)
# print(batata2.uk)

nm.solve_normal_system(batata2)
print("Input:\t[2.3 \t\t3.7 \t\t0.3 \t4.2]")
print("solution:{}".format(batata2.a))

error = nm.quatratic_error(batata2)
print("error:{}".format(error))
