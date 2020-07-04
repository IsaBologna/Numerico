import numpy as np
import sys
import math

from ep_methods.tarefas import Teste, Teste_B, Teste_C
import ep_methods.numerical_methods as nm

# N=128
# T=1
# p=np.array([0.35])

# batata = Teste(N,T,p)

# nm.matrix_uk(batata)
# batata.exact_solution() #calcula uT a partir de uk
# nm.solve_normal_system(batata)
# print(batata.a)

# p2=np.array([0.15,0.3,0.7,0.8])

# batata2 = Teste_B(N,T,p2)

# nm.matrix_uk(batata2)
# batata2.exact_solution()
# nm.solve_normal_system(batata2)
# print("Input:\t[2.3 \t\t3.7 \t\t0.3 \t4.2]")
# print("solution:{}".format(batata2.a))

# error = nm.quatratic_error(batata2)
# print("error:{}".format(error))


#* Leitura dos dados do arquivo .txt
uT = []
with open('test.txt') as f:
    first_line = f.readline()
    l1=first_line.split('       ') 
    
    p=np.zeros(len(l1))
    for i in range(0,len(l1)):
        p[i] = float(l1[i]) # vetor com os pontos para a fonte de calor
    
    for line in f: 
        uT.append(float(line)) # uT(xi)

N=128
T=1

batata=Teste_C(N,T,p)
nm.matrix_uk(batata)
batata.exact_solution(uT) #calcula uT a partir de uk
nm.solve_normal_system(batata)
print("solution:{}".format(batata.a))

error = nm.quatratic_error(batata)
print("error:{}".format(error))