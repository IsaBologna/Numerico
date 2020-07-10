#!/usr/bin/env python3

import numpy as np
import sys
import math
from prompt_toolkit import prompt

from ep_methods.tarefas import Teste_A, Teste_B, Teste_C, Teste_D
import ep_methods.numerical_methods as nm


if __name__ == '__main__':
    T = 1
    try:
        while True:
            choice = ''
            item = None
            while choice not in ['a', 'b', 'c', 'd']:
                choice = prompt("\n\nDeseja executar qual item? [a, b, c, d]: ")

                if choice == 'a':
                        N = 128
                        p = np.array([0.35])
                        item = Teste_A(N, T, p)

                        nm.matrix_uk(item)
                        item.exact_solution()
                        nm.solve_normal_system(item)
                        print("Intensidades ak: {}".format(item.a))

                        error = nm.quatratic_error(item)
                        print("Erro: {}".format(error))
                        
                elif choice == 'b':
                    N = 128
                    p = np.array([0.15, 0.3, 0.7, 0.8])
                    item = Teste_B(N, T, p)

                    nm.matrix_uk(item)
                    item.exact_solution()
                    nm.solve_normal_system(item)
                    print("Intensidades ak: {}".format(item.a))

                    error = nm.quatratic_error(item)
                    print("Erro: {}".format(error))

                elif choice == 'c':
                    #* Leitura dos dados do arquivo .txt
                    uT = []
                    with open('test.txt') as f:
                        first_line = f.readline()
                        l1 = first_line.split('       ') 
                
                        p = np.zeros(len(l1))
                        for i in range(0,len(l1)):
                            p[i] = float(l1[i]) # vetor com os pontos para a fonte de calor
                        
                        for line in f: 
                            uT.append(float(line)) # uT(xi)

                    N = int(prompt("Digite o valor de N: "))
                    item = Teste_C(N, T, p)
                    
                    nm.matrix_uk(item)
                    item.exact_solution(uT)
                    nm.solve_normal_system(item)
                    print("Intensidades ak: {}".format(item.a))

                    error = nm.quatratic_error(item)
                    print("Erro: {}".format(error))

                    nm.plot_solution(item)

                elif choice == 'd':
                    #* Leitura dos dados do arquivo .txt
                    uT = []
                    with open('test.txt') as f:
                        first_line = f.readline()
                        l1=first_line.split('       ') 
                
                        p=np.zeros(len(l1))
                        for i in range(0,len(l1)):
                            p[i] = float(l1[i]) # vetor com os pontos para a fonte de calor
                        
                        for line in f:
                            ruido = 1 + 0.01 * ((np.random.random() - 0.5) * 2)  # Ruido = 1 + 0.01r, r ∈ [-1, 1]
                            uT.append(ruido * float(line))  # uT(xi) * ruido

                    N = int(prompt("Digite o valor de N: "))
                    item = Teste_D(N, T, p)
                    
                    nm.matrix_uk(item)
                    item.exact_solution(uT)
                    nm.solve_normal_system(item)
                    print("Intensidades ak: {}".format(item.a))

                    error = nm.quatratic_error(item)
                    print("Erro: {}".format(error))
                    
                    nm.plot_solution(item)

                else:
                    print("Opcao invalida :(")
    except (KeyboardInterrupt, EOFError):
        print("Bye bye")
        sys.exit()
