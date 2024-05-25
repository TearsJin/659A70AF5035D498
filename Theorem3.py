import time
from Crypto.Util.number import getPrime,isPrime
from tqdm import tqdm
from sage.all import *

load("Lemma3.sage")

def binchange(num, base):
    s = ''
    while num:
        i = num % base
        s += str(i)
        num = int(num / base)
    return s[::-1]

def lift(f, df, p, k, previous):
    result = []
    for lower_solution in previous:
        dfr = df(lower_solution)
        fr = f(lower_solution)
        if dfr % p != 0:
            t = (-(xgcd(dfr, p)[1]) * (fr // p ** (k - 1))) % p
            result.append(lower_solution + t * p ** (k - 1))
        if dfr % p == 0:
            if fr % p ** k == 0:
                for t in range(0, p):
                    result.append(lower_solution + t * p ** (k - 1))
    return result

def hensel_lifting(f, df, p, k, base_solution):
    if type(base_solution) is list:
        solution = base_solution
    else:
        solution = [base_solution]
    for i in range(2, k + 1):
        solution = lift(f, df, p, i, solution)
    return solution

def Generate(k,lb,Bq,r,beta,gamma):
    Qbits = (k // 2 * lb + Bq + 2)
    Nbits = int(Qbits / beta)
    Pbits = Nbits - Qbits * r
    Nbits = Qbits * r + Pbits

    Bp = Pbits - (k // 2 * lb + 2)

    print("[+] Bq: {}, Bp: {}".format(Bq,Bp))
    while True:
        Sq = [getPrime(lb - 1) for _ in range(k // 2 - 1)]
        a = getPrime(Bq - 1)
        i = 0
        while i < 100:
            q_ = getPrime(lb)
            Q = int(2 * a * prod(Sq) * q_ + 1)
            if isPrime(Q):
                break
            i += 1
        if i == 100:
            continue
        
        Sp = [getPrime(lb - 1) for _ in range(k // 2 - 1)]
        b = getPrime(Bp - 1)
        
        i = 0
        while i < 100:
            p_ = getPrime(lb)
            P = int(2 * b * prod(Sp) * p_ + 1)
            if isPrime(P):
                break
            i += 1
        if i == 100:
            continue
        
        N = P * Q ** r
        E1 = [q_]
        E2 = [p_]
        while RR(log(prod(E1) * prod(E2),N)) <= gamma:
            p = Sp.pop(0)
            if p not in E2:
                E1.append(p)
            
            q = Sq.pop(0)
            if q not in E1:
                E2.append(q)
        E1 = E1[1:]
        E2 = E2[1:]

        assert gcd(P - 1,Q - 1) == 2
        assert (P - 1) % prod(E1) == 0 
        assert (Q - 1) % prod(E2) == 0 
        return (P,Q,prod(E1),prod(E2))


def Theorem3(N,E1,E2,r,beta = None):
    e = E1 * E2
    gamma = RR(log(e,N))

    GCD = gcd(E1,E2)
    s1 = (inverse_mod(E1 // GCD,E2 // GCD) * ((N - 1) // GCD)) % E2
    s = E1 * s1 + 1
    ur = Integer((N * inverse_mod(s,e)) % e)
    es = factor(e)
    print(f'[+] es = {es}')
    print("[+] ur =",ur)
    candidateU = []
    eiz = []       # e_i
    for ei in es:
        PR = PolynomialRing(Zmod(ei[0]),"x")
        x = PR.gens()[0]
        f = x ** r - (ur % ei[0])
        print("[+] ur1:",(ur % ei[0]))
        U = [int(root[0]) for root in f.roots()]
        if ei[1] == 1:
            candidateU.append(U)
        else:
        # Hensel lifting algorithm
            PR = PolynomialRing(ZZ,"x")
            x = PR.gens()[0]
            f = x ** r - (ur % (ei[0] ** ei[1]))
            print("[+] ur2:",(ur ,ei[0] , ei[1]))
            candidateU.append(hensel_lifting(f,derivative(f,x),ei[0],ei[1],U))
            for _ in candidateU[-1]:
                print(_,f)
                print(f(_) % ei[0] ** ei[1])
        eiz.append(int(ei[0] ** ei[1]))
    assert prod(eiz) == e

    Us = []
    for i in range(r ** len(es)):
        try:
            Us.append([U[int(j)] for U,j in zip(candidateU,binchange(i,r).zfill(len(es)))])
        except:
            pass
    print(f'[+] The number of candidate u: {len(Us)}')
    R = Integers(N)
    PR = PolynomialRing(R,"x",1)
    x = PR.gens()[0]
    for us in tqdm(Us): 
        u = CRT(us,eiz) # Chinese Remainder Theorem
        f = x + (u * inverse_mod(e,N))
        if beta == None:
            mu = gamma - (1 / (4 * r)) - 0.01
            print("[+] mu:",mu)
            print(f'[+] gamma - mu > 1 / 4r : {RR(gamma - mu) >= RR(1 / 4 / r)}')
            result = iterationAlgorithm(f,r,gamma,mu)
            if result != []:
                return list(set([gcd(r + (u * inverse_mod(e,N)) ,N) for r in result]))
        else:
            mu = gamma - (beta - r * beta ** 2) - 0.01
            assert RR(gamma) >= RR((beta - r * beta ** 2))
            result = Coppersmith(f,r,1,beta,mu)
            if result != []:
                return gcd(result[0][0] + (u * inverse_mod(e,N)),N)

gamma = 0.16                     # e = N ^ gamma
beta = 0.2                       # Q = N ^ beta
r = 3                            # N = P * Q ^ r
k = 16                           # The number of prime
lb = 12                          # bits-length of primes
Bq = 50                          # bits-length of large prime in P

RR = RealField(2 ** 5)
P,Q,E1,E2 = Generate(k,lb,Bq,r,beta,gamma)
N = P * Q ** r

gamma = RR(Integer(E1 * E2).nbits() / Integer(N).nbits())
beta = RR(Integer(Q).nbits() / Integer(N).nbits())

print("[+] beta: {}, gamma: {}".format(beta,gamma))
Starttime = time.time()
# print(Theorem3(N,E1,E2,r))
print(Theorem3(N,E1,E2,r,beta))
print("[+] runtime:",time.time() - Starttime)