#!/usr/bin/env python3

import sys
from prompt_toolkit import prompt

from ep1_methods import Item_A1, Item_A2, Item_B, Item_C
import ep1_methods.numerical_methods as nm


T = 1
try:
    while True:
        tarefa = prompt("Executar primeira ou segunda tarefa? [1/2]: ")

        if tarefa == '1':
            N = int(prompt("Digite o valor de N: "))
            Lambda = float(prompt("Digite o valor de lambda: "))

            M = int((T * N**2) / Lambda)

            choice = ''
            while choice not in ['a', 'b', 'c']:
                choice = prompt("Deseja executar qual item? [a, b, c]: ")

                if choice == 'a':
                    item = Item_A2(N, M, Lambda, T)
                elif choice == 'b':
                    item = Item_B(N, M, Lambda, T)
                elif choice == 'c':
                    item = Item_C(N, M, Lambda, T)
                else:
                    print("Opcao invalida :(")

            nm.finite_difference_method(item)
            nm.plot(item.u, item.dx, item.M, title="Solução item {}".format(choice.upper()))

            while True:
                cmd = prompt("\nPressione enter para voltar ao inicio... ")
                break


        elif tarefa == '2':
            N = int(prompt("Digite o valor de N: "))
            Lambda = N
            M = N

            choice = ''
            while choice not in ['a', 'b', 'c']:
                choice = prompt("Deseja executar qual item? [a, b, c]: ")

                if choice == 'a':
                    item = Item_A2(N, M, Lambda, T)
                elif choice == 'b':
                    item = Item_B(N, M, Lambda, T)
                elif choice == 'c':
                    item = Item_C(N, M, Lambda, T)
                else:
                    print("Opcao invalida :(")
            
            cmd = prompt("Resolver com Euler ou Crank-Nicolson? [E / C]: ")
            if cmd == 'E':
                A, B = nm.create_A_matrix(Lambda, N-1)
                L, D = nm.LDL_decomposition(A, B)

                nm.implicit_euler_method(item, L, D)
                nm.plot(item.u, item.dx, item.M, title="Solução item {}".format(choice.upper()))
            elif cmd == 'C':
                A, B = nm.create_A_matrix(Lambda/2, N-1)
                L, D = nm.LDL_decomposition(A, B)

                nm.crank_nicolson_method(item, L, D)
                nm.plot(item.u, item.dx, item.M, title="Solução item {}".format(choice.upper()))
            else:
                print("Opcao invalida :(")

            
            while True:
                cmd = prompt("\nPressione enter para voltar ao inicio... ")
                break


        else:
            print("Opcao invalida :(")

except (KeyboardInterrupt, EOFError):
    print("Bye bye")
    sys.exit()
