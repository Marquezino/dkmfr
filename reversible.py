"""reversible.py: Description ..."""

__author__ = "E.Dalcumune, L.Kowada, F.Marquezino, C.Figueiredo, A.Ribeiro"
__email__ = "franklin@cos.ufrj.br"
__licence__ = "MIT"


from math import log
from sympy.combinatorics.permutations import Permutation


Ltfc = []
LRtfc = []


def matrix_perm(n):
    """ The function generates a matrix of dimension 3^(n-1)*n and stores its
    elements to be use later.

    Args:
        n (int): number of bits.

    Returns:
        matrix: matrix with all possible logic gates for the Generalized
        Toffoli library for n-bit circuits.
        n_res: combinations for all lines."""

    seq = [0, 2, 3]                 # number of possibilities for each line
    # combinations for non target lines
    res = [[0 for j in range(n-1)] for i in range(3**(n-1))]
    # combinations for all lines
    n_res = [[0 for j in range(n-1)] for i in range(3**(n-1))]
    pos = [0]*(3**(n-1))            # change in the number of controls

    for i in range(3**(n-1)):
        k = i
        for j in range(n-1):
            res[i][j] = seq[k % 3]
            k = k//3

    position = [0]*(n)
    for c in range(1, n):           # c is the number of controls
        d = (combination(n-1, c-1) * (2**(c-1)))
        # change position in the number of controls
        position[c] = (position[c - 1] + int(d))

    ctrl = [0]*3**(n-1)
    # new position according to the number of controls
    for i in range(3**(n-1)):
        c = controles(res[i], n-1)
        ctrl[i] = c
        pos[i] = position[c]
        position[c] += 1

    for i in range(3**(n-1)):
        n_res[pos[i]] = res[i]

    matrix = []
    for i in range(3**(n-1)):
        matrix.append([0]*n)
    # lines for totally controlled gates n-1 controls
    for i in range(3**(n-1)-1, 3**(n-1)-2**(n-1)-1, -1):
        for k in range(0, n, 1):
            seq_temp = n_res[i][:]
            seq_temp.insert(k, 1)
            matrix[i][k] = seq2perm(seq_temp, n)
    # lines for gates with n-2 controls until 0 control
    for i in range(3**(n-1) - 2**(n-1)-1, -1, -1):
        pos_m = 0
        pos_zero = -1
        for j in range(n-1):
            pos_m += ((n_res[i][j])+1)//2*(3**j)
            if n_res[i][j] == 0:
                pos_zero = j
        pos_p = pos_m + 3**pos_zero
        pos_q = pos_m + 2*3**pos_zero
        p = pos[pos_p]
        q = pos[pos_q]
        for j in range(n):
            matrix[i][j] = matrix[p][j] * matrix[q][j]

    return matrix, n_res


def strGates(n, a, b):
    """.

    Args:
        n:
        a:
        b:

    Returns:
        texto: ."""

    texto = 'T' + str(n) + ' '
    c = a ^ b
    for i in range(n):
        if c & 1:
            target = i
        else:
            if (a >> i) & 1:
                texto += "b" + str(i) + ","
            else:
                texto += "b" + str(i) + "',"
        c = c >> 1
    texto += "b" + str(target)
    return texto


def calc_dist(pi, countGates):
    """ The function manages a set of functions that are designed to find
    suitable gates for the circuit.

    Args:
        pi: the permutation resulting from the last application of gates.
        countGates: number of gates already found.

    Returns:
        coutGates: updated number of gates already found.
        gates: gates applied."""

    n = int(log(pi.size, 2))
    size = pi.size
    iota = Permutation([], size=2**n)
    gates = []
    portas = []

    while pi != iota:
        atempt, porta = search_2move(pi)
        if atempt is not None:
            pi = atempt
            gates.append(porta)
            countGates += 1
            continue

        atemp, aux, portas = search_seq(pi)
        if aux != 0:
            pi = atemp
            gates += portas
            countGates += aux
            continue

        atempt, porta = search_0move(pi)
        if atempt is not None:
            pi = atempt
            gates.append(porta)
            countGates += 1
            continue

        atempt, count, portas = Alg_2(pi)
        if atempt is not None:
            pi = atempt
            gates += portas
            countGates += count
            continue


        print('Number of gates:', ' ====> ', countGates)
        print('***************** Could not solve ****************', pi, '\n')
        return 0
    return countGates, gates


def search_2move(pi):
    ''' The function tries to find a gate that when applied to the permutation
    finds a new permutation such that
    DH(permutation,identity)=DH(new permutation,identity)+2.

    Args:
        pi: the permutation resulting from the last application of gates.

    Returns:
        pi_return: new permutation.
        T_temp: gate found.'''

    iota = Permutation([], size=pi.size)
    T_temp = iota
    n = int(log(pi.size, 2))
    dmin = pi.size
    pi_return = None

    chosen_element = None
    chosen_neighbor = None

    for cycle in pi.cyclic_form:
        for element in cycle:
            for i in range(len(cycle)):
                if element == cycle[i]:
                    d_atual = dh_int_int(cycle[i-1], cycle[i])
            for neighbor in neighbors(element, n):
                piaux = pi*swap(element, neighbor)
                # checks if there is a 2-move
                if dh_perm_perm(piaux) < dh_perm_perm(pi):
                    # checks if it is the 1st 2-move found
                    if pi_return is None:
                        pi_return = piaux
                        chosen_element = element
                        chosen_neighbor = neighbor
                    # checks if the 2-move joins cycles and if the H.D.
                    # between treated neighbors is minimal
                    if P(piaux) > P(pi) and d_atual < dmin:
                        dmin = d_atual
                        pi_return = piaux
                        chosen_element = element
                        chosen_neighbor = neighbor

    if chosen_element is not None and chosen_neighbor is not None:
        # saves gate to ".tfc" file
        Ltfc.append(strGates(n, chosen_element, chosen_neighbor))
        # returns gate to write in ".txt" file
        T_temp = swap(chosen_element, chosen_neighbor)

    return pi_return, T_temp


def search_seq(pi):
    ''' The function checks whether, in a cycle, there is a sequence of
    0-moves that ends with a 2-move and applies the sequence of gates.

    Args:
        pi: the permutation resulting from the last application of gates.

    Returns:
        pi_return: new permutation.
        countGates: updated number of gates already found.
        gates: gates in the sequence.'''

    n = int(log(pi.size, 2))
    for cycle in pi.cyclic_form:
        gates = []
        countGates = 0
        flagServe = 0
        if len(cycle) == 2:  # if there was 2-move, it had been applied before
            break
        beginSeq = len(cycle)       # mark first distance > 1
        endSeq = len(cycle)         # mark second distance > 1
        for i in range(len(cycle)):
            if dh_int_int(cycle[i-1], cycle[i]) != 1:
                if beginSeq == len(cycle):  # first time in the cycle dist > 1
                    beginSeq = i-1
                elif endSeq == len(cycle):  # second time in the cycle dist > 1
                    endSeq = i-1
                if endSeq < len(cycle):    # there are at least 2 with dist > 1
                    # to ensure that the last move is 2-move
                    if dh_int_int(cycle[beginSeq], cycle[endSeq]) == 1:
                        piaux = pi
                        for j in range(beginSeq+1, endSeq):
                            # not necessarily 0-move
                            piaux = piaux*swap(cycle[j], cycle[j+1])
                            countGates += 1
                            # returns gate to write to ".txt" file
                            gates.append(swap(cycle[j], cycle[j+1]))
                            Ltfc.append(strGates(n, cycle[j], cycle[j+1]))
                        # does not make a possible move
                        return piaux, countGates, gates
                    else:
                        beginSeq = endSeq        # start looking again
                        # do not use endSeq = len (cycle) to reset it
                        endSeq = len(cycle)+1
        # if all distances between neighbors are 1,
        # then disassemble the entire cycle
        if beginSeq == len(cycle):
            piaux = pi
            for i in range(len(cycle)-1):
                piaux = piaux*swap(cycle[i-1], cycle[i])
                countGates += 1
                # returns gate to write to ".txt" file
                gates.append(swap(cycle[i-1], cycle[i]))
                Ltfc.append(strGates(n, cycle[i-1], cycle[i]))
            return piaux, countGates, gates
        # if there is only a distance other than 1,
        # then disassemble the entire cycle
        elif endSeq == len(cycle):
            piaux = pi
            for i in range(beginSeq+1, len(cycle)-1):
                piaux = piaux*swap(cycle[i], cycle[i+1])
                countGates += 1
                # returns gate to write to ".txt" file
                gates.append(swap(cycle[i], cycle[i+1]))
                Ltfc.append(strGates(n, cycle[i], cycle[i+1]))
            for i in range(beginSeq+1):
                piaux = piaux*swap(cycle[i-1], cycle[i])
                countGates += 1
                # returns gate to write to ".txt" file
                gates.append(swap(cycle[i-1], cycle[i]))
                Ltfc.append(strGates(n, cycle[i-1], cycle[i]))
            return piaux, countGates, gates
    return None, 0, None


def search_0move(pi):
    ''' The function tries to find a 0-move that joins cycles.

    Args:
        pi: the permutation resulting from the last application of gates.

    Returns:
        pi_return: new permutation.
        gate: 0-move found.'''

    n = int(log(pi.size, 2))
    dmin = pi.size
    pi_return = None
    gate = None

    chosen_element = None
    chosen_neighbor = None

    for cycle in pi.cyclic_form:
        for element in cycle:
            for i in range(len(cycle)):
                if element == cycle[i]:
                    d_atual = dh_int_int(cycle[i-1], cycle[i])
            for neighbor in neighbors(element, n):
                piaux = pi*swap(element, neighbor)
                if dh_perm_perm(piaux) == dh_perm_perm(pi) and P(piaux)< P(pi):
                    if d_atual < dmin:
                        dmin = d_atual
                        pi_return = piaux
                        # gate to write to ".txt" file
                        gate = swap(element, neighbor)
                        chosen_element = element
                        chosen_neighbor = neighbor
    # saves gate to ".tfc" file
    if chosen_element is not None and chosen_neighbor is not None:
        Ltfc.append(strGates(n, chosen_element, chosen_neighbor))

    return pi_return, gate


def Alg_2(pi):
    ''' The function checks whether there is a sequence of 0-moves that
    ends with a 2-move and applies the sequence of gates.

    Args:
        pi: the permutation resulting from the last application of gates.

    Returns:
        piaux: new permutation.
        count: updated number of gates already found.
        gates: gates found.'''

    size = pi.size
    piaux = []
    gates = []
    count = 0
    flag = 0
    j = size - 1
    while j > 0 and flag == 0:
        while j == pi(j):
            j -= 1
        piaux, count, flag, gates = replace(pi, j, pi(j), count, gates)
        pi = piaux
        j -= 1
    return piaux, count, gates


def replace(pi, j, i, count, gates):
    '''The function updates i and i_line while the Hamming distance between
    i and j is different from 1, then applies the gate and returns to the Alg_2
    function.

    Args:
        pi: the permutation resulting from the last application of gates.
        j,i: elements of the permutation.
        count: number of gates already found.
        gates: gates already found.

    Returns:
        pi_return: new permutation.
        count: updated number of gates already found.
        dm: flag.
        gates: gates found.'''

    n = int(log(pi.size, 2))
    indice = 0
    pi_return = []
    i_line = 0
    # While it is not possible to put j in its position, do (i = i')
    while dh_int_int(i, j) != 1:
        x = i ^ j
        k = 0
        # signal is equal to 1 whenever bit i_k is different from j_k
        while k < n:
            zeroum = (x >> k) & 1
            if zeroum == 1:
                indice = k
                break
            k += 1
        i_line = i ^ (2**(indice))
        piaux = pi * swap(i, i_line)
        count += 1
        gates.append(swap(i, i_line))  # returns gate to write to ".txt" file
        Ltfc.append(strGates(n, i, i_line))
        i = i_line
        pi = piaux
    pi_return = pi*swap(i, j)
    count += 1
    gates.append(swap(i, j))   # returns gate to write to ".txt" file
    Ltfc.append(strGates(n, i, j))
    if dh_perm_perm(pi_return) < dh_perm_perm(pi):
        dm = 1
    else:
        dm = 0
    return pi_return, count, dm, gates


def inv(n, perm):
    ''' The function finds the inverse permutation.

    Args:
        n: number of bits.
        perm: a permutation.

    Returns:
        inverse: the inverse permutation.'''

    inverse = [0]*2**n
    for i, p in enumerate(perm):
        inverse[p] = i
    return inverse


def seq2perm(seq, n):
    ''' The function turns a sequence into a permutation.

    Args:
        seq: a sequence.
        n: number of bits.

    Returns:
        p: resulting permutation.'''

    a = b = 0
    p = Permutation([], size=2**n)
    for i in range(n):
        if seq[i] == 1:
            a += 0*2**i
            b += 1*2**i
        elif seq[i] == 2:
            a += 0*2**i
            b += 0*2**i
        elif seq[i] == 3:
            a += 1*2**i
            b += 1*2**i
    p = Permutation(a, b)
    return p


def controles(x, n):
    ''' The function returns the number of controls for the gate.

    Args:
        x: ternary number.
        n: number of bits.

    Returns:
        c: number of controls.'''

    c = 0
    for i in range(n):
        if x[i]:   # checks if it is different from zero
            c += 1
    return c


def fatorial(n, limite=1):
    ''' The function computes the factorial number.

    Args:
        n: number of bits.

    Returns:
        resultado: factorial number.'''

    if n == 0 or n == 1:
        return 1
    if limite < 1 or limite > n:
        return -1
    else:
        fatArr = range(limite, n + 1)
        resultado = 1
        for i in fatArr:
            resultado = resultado * i
        return resultado


def combination(n, p):
    ''' The function computes the combination C^{n}_{p}.

    
    Args:
        n,p:

    Returns:
        a/b: the combination C^{n}_{p}.'''

    if p > n:
        return -1
    a = b = limite = 1
    np = n - p
    if np < p:
        limite = np+1
        b = p
    else:
        limite = p+1
        b = np
    a = fatorial(n, limite)
    b = fatorial(b)
    return a/b


def P(perm):
    ''' The function computes the P() function described in the paper, which
    considers the S() function and the size of the cycle. The result of
    the computation serves as a parameter for choosing the gate using the
    function 2-move().

    Args:
        perm: the permutation resulting from the last application of gates.

    Returns:
        aux: parameter found.'''

    aux = 0
    for cycle in perm.cyclic_form:
        aux += S(cycle)/len(cycle)
    return aux


def S(cycle):
    ''' The function computes the sum of Hamming distances between
    neighboring elements of a cycle.

    Args:
        cycle: cycle in a permutation.

    Returns:
        aux: sum of Hamming distances.'''

    size = len(cycle)
    aux = 0
    for i in range(size):
        aux += dh_int_int(cycle[i-1], cycle[i])
    return aux


def neighbors(a, size):
    ''' The function discovers the elements of a permutation that are neighbors
    (DH = 1) of a given element.

    Args:
        cycle: cycle in a permutation.

    Returns:
        aux: sum of Hamming distances.'''

    aux = []
    for i in range(size):
        aux.append(a ^ (1 << i))
    return aux


def dh_perm_perm(p, q=None):
    ''' The function computes the Hamming distance between two given
    permutations.

    Args:
        p,q: two permutations.

    Returns:
        dist: Hamming distance between p and q.'''

    if q is None:
        q = Permutation([], size=p.size)
    elif p.size != q.size:
        raise ValueError('permutations must be of same sizes')

    dist = 0
    for i in range(p.size):
        pi = p.array_form[i]
        qi = q.array_form[i]
        dist += dh_int_int(pi, qi)

    return dist


def dh_perm_int(p, a):
    ''' The function computes the Hamming distance between an element
    of the permutation p and the integer a.

    Args:
        p: permutation.
        a: integer.

    Returns:
        Hamming distance between p.array_form[a] and a.'''

    return dh_int_int(p.array_form[a], a)


def dh_int_int(a, b, bits=64):
    ''' The function computes the Hamming distance between two integer
    numbers.

    Args:
        a,b: integer numbers.

    Returns:
        Hamming distance between a and b.'''

    x = a ^ b
    return sum((x >> i & 1) for i in range(bits))


def swap(i, j):
    ''' The function makes the swap between i and j, if the Hamming
    distance between them is equal to 1.

    Args:
        i,j: elements in a permutation.

    Returns:
        swap between elements i and j.'''

    if dh_int_int(i, j) == 1:
        return Permutation(i, j)
    else:
        raise ValueError('Invalid Transposition!!')
        return 0


def unicyclic(n, index=0):
    ''' .

    Args:
        n: number of bits.
        index:

    Returns:
        .'''

    sigma = Permutation([], size=(2**n) - 1)
    sigma = sigma + index
    return Permutation(list(Permutation([sigma.list(2**n)])))
