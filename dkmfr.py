#!/usr/bin/env python

""" dkmfr.py: Experiments with reversible circuit synthesis."""

__author__ = "E.Dalcumune, L.Kowada, F.Marquezino, C.Figueiredo, A.Ribeiro"
__email__ = "franklin@cos.ufrj.br"
__licence__ = "MIT"


import sys
import random
import copy
import time
import benchmarks as bm
import reversible as rev
import argparse
import itertools


from math import log
from math import factorial
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.permutations import Cycle

from pathlib import Path



def str_time():
    return time.strftime("%Y-%m-%d_%H-%M-%S",time.gmtime())


def main(args=None):

    str_current_time = str_time() #string corresponding to the current time
    Path("output/" + str_current_time).mkdir(parents=True, exist_ok=True)  #create folder output if it does not exist
    out_path = 'output/' + str_current_time + '/'

    f_time = open(out_path + 'Time.txt', mode = 'w') #open file to save the elapsed time for each permutation

    gs = open(out_path + 'Number_Gates_', mode='w')       # abre arquivo p/ imprimir distâncias p/ pi  
    
    start_time = time.time()

    parser = argparse.ArgumentParser(description='''example: $ DKMFR.py
                                     --bi 2 --functions hwb4 4b15g_2''')
    parser.add_argument('--bi', action='store', dest='bi', default=2,
                        required=False, help='choose 1 for forward only, and '
                        ' choose 2 for forward/reverse.', type=int)
    parser.add_argument('--functions', nargs='+', action='store',
                        dest='functions', help='''Choose benchmark functions
                        from 4b15g_1, 4b15g_2, 4b15g_3, 4b15g_4, 4b15g_5,
                        hwb4, 4_49, toffoli_double_2, nth_prime4_inc,
                        ham3_complete_47(28), miller_complete_5, toffoli_1,
                        ex-1_82, 3_17, nth_prime3_inc, RCFK, ex1Miller,
                        ex2Miller, ex3Miller, ex4Miller, ex5Miller,
                        ex6Miller, ex7Miller, aj-e11_complete_74(81),
                        mod5mils_complete_26(18), nth_prime5_inc, hwb5_13,
                        mod5adder, hwb6, nth_prime6_inc,
                        graycode6_complete_19, nth_prime7_inc, ham7,
                        hwb7_15, random. You can choose multiple functions
                        separating them by spaces.''',
                        required=False, default=['hwb4'])
    parser.add_argument('--samples', action='store', dest='samples',
                        default=10, type=int,
                        help='how many random permutations will be generated.')
    parser.add_argument('--nbits', action='store', dest='nbits',
                        default=3, type=int,
                        help='how many bits in random permutations.')
    args = parser.parse_args(args)

    func_list = [args.samples*[f] if f == 'random' else [f] for f in args.functions]
    func_list = list(itertools.chain(*func_list))
    print(func_list)

    random_seq = 0  # this counter will be used to name the output files

    for option in func_list:
        if option == 'random':
            temp = list(range(2**args.nbits))
            random.shuffle(temp)
            pi = Permutation(temp)
            random_seq += 1
        else:
            pi = bm.pi_dict[option]

        # files to print gates for pi
        if option == 'random':
            file_name = 'S' + str(args.bi) + '-' + option + '-' + str(random_seq)
        else:
            file_name = 'S' + str(args.bi) + '-' + option
        arc = open(out_path + file_name+'.txt', mode='w')
        arc_tfc = open(out_path + file_name+'.tfc', mode='w')

        size = pi.size                        # size of permutation
        iota = Permutation([], size=pi.size)  # identity permutation
        n = int(log(size, 2))                 # number of bits

        info = '# permutation: ' + str(option)
        info += '\n# ' + str(pi.array_form)
        info += '\n# n.bits: ' + str(n)
        info += '\n# DH(pi)= ' + str(rev.dh_perm_perm(pi))
        info += '\n# reverse option: ' + str(args.bi-1)
        print(info, file=arc)
        print(info, file=arc_tfc)

        # header for .tfc gates files
        c = ['.v ', '\n.i ', '\n.o ']
        for j in range(3):
            arc_tfc.write(c[j])
            for k in range(n-1):
                print("b%d" % k, end=',', file=arc_tfc)
            print("b%d" % (n-1), end=' ', file=arc_tfc)

        print('\nBEGIN', file=arc_tfc)

        print('----------------------------------------------')
        print(option)
        print('DH(pi) = ', rev.dh_perm_perm(pi))

        # Matrix obtained in pre-processing
        matrix, lines = rev.matrix_perm(n)
        end_time = time.time()

        # 2 - One permutation or permutations in sequence

        # pi = iota + 1
        pi_final = iota

        start_time = time.time()         # time for pi in seconds

        while pi != pi_final:       # for permutations in sequence
            fi = pi
            countGates = 0          # number of Gates
            nc = 0                  # number of controls
            start = 0               # matrix line corresponding to change of nc
            LD = []
            LR = []
            LRtfc = []

            while nc < (n-1):
                for reverse in range(args.bi):
                    temp = fi
                    signal = 0
                    # number of possible gates with nc controls
                    lim = (rev.combination(n-1, nc) * (2**(nc)))
                    for i in range(int(lim)):
                        # number of matrix columns
                        for j in range(n):
                            if reverse:
                                piaux = matrix[start + i][j] * fi
                            else:
                                piaux = fi * matrix[start + i][j]
                            if rev.dh_perm_perm(piaux)< rev.dh_perm_perm(temp):
                                temp = piaux
                                signal = 1
                                gate = matrix[start+i][j]
                                controls = start+i
                                target = j
                    if signal == 0:  # if no more gates with nc controls
                        nc += 1
                        start += int(lim)
                    else:
                        countGates += 1  # if more gates with nc controls
                        if reverse:
                            LD.append(gate)
                        else:
                            LR.append(gate)

                        # for write gates to .tfc file
                        toffoli = 1
                        for k in range(n-1):
                            if lines[controls][k] != 0:
                                toffoli += 1
                        line_tfc = 'T' + str(toffoli) + ' '

                        # negative / positive controls
                        for k in range(n):
                            shift_flag = k > target
                            if k == target:
                                continue
                            if lines[controls][k-shift_flag] == 2:
                                line_tfc += "b" + str(k) + "',"
                            if lines[controls][k-shift_flag] == 3:
                                line_tfc += 'b' + str(k) + ","

                        line_tfc += 'b' + str(target)
                        if reverse:
                            rev.LRtfc.append(line_tfc)
                        else:
                            rev.Ltfc.append(line_tfc)

                    fi = temp

            gatesPi = []
            dist, gatesPi = rev.calc_dist(fi, countGates)

            if dist == 0 and pi != iota:
                print('Error!! Unable to solve ', file=f, flush=True)
                print('Error!! Unable to solve ', pi)
                break

            print('Total number of gates = ', dist)
            print(dist, file=gs) # Adicionado para resolver 2
            for gtotal in gatesPi:
                LR.append(gtotal)
                
            print('elapsed time = ', time.time() - start_time, file = f_time)
        
        #    pi = pi+1
            pi = pi_final

        for lineL in rev.LRtfc:
            print(lineL, file=arc_tfc)
        for lineL in rev.Ltfc[::-1]:
            print(lineL, file=arc_tfc)

        print('END', file=arc_tfc)
        end_time = time.time()
        
        
#        print('elapsed time = ', end_time - start_time)

       
        print('# Total Gates = ', dist, file=arc)
        print('# Total Gates = ', dist, file=arc_tfc)

        for x in LD[::-1]:
            LR.append(x)
        print('Gates Sequence:', file=arc)
        for x in LR:
            print(x, file=arc)

        rev.LRtfc = [] # Adicionado para resolver 3 
        rev.Ltfc = [] # Adicionado para resolver 3
            


if __name__ == "__main__":
    sys.exit(main())
