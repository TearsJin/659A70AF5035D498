import time
from Crypto.Util.number import getPrime,isPrime
from sage.all import *

load("Lemma3.sage")

def Generate(k,lb,Bq,r,beta,gamma):
    Qbits = (k // 2 * lb + Bq + 2)
    Nbits = int(Qbits / beta)
    Pbits = Nbits - Qbits * r
    Nbits = Qbits * r + Pbits

    print(Nbits * gamma)
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
        e = getPrime(int(Nbits * gamma) - 1) # log(e * b,N) = Bp 
        b = getPrime((Bp - int(Nbits * gamma)) - 1)
        
        i = 0
        while i < 100:
            p_ = getPrime(lb)
            P = int(2 * b * e * prod(Sp) * p_ + 1)
            if isPrime(P):
                break
            i += 1
        if i == 100:
            continue

        print(log(P,2).n(),log(Q,2).n())
        return (P,Q,e)

def Theorem1(N,e,r,beta = None):
    gamma = RR(log(e,N))
    PR = PolynomialRing(Zmod(e),"x")
    x = PR.gens()[0]
    f = x ** r - (N % e) 
    TTT = time.time()
    us = [int(root[0]) for root in f.roots()]
    print("[+] AMM time :", time.time() - TTT)

    PR = PolynomialRing(Zmod(N),"x",1)
    x = PR.gens()[0]
    print("[+] len of us",len(us))
    for u in us:
        f = x + (u * inverse_mod(e,N))
        if beta == None:
            # beta is unknown
            mu = gamma - (1 / (4 * r)) - 0.01
            print("[+] mu:",mu)
            result = iterationAlgorithm(f,r,gamma,mu)
            if result != []:
                return list(set([gcd(r + (u * inverse_mod(e,N)) ,N) for r in result]))
        else:
            # beta is known
            mu = gamma - (beta - r * beta ** 2) - 0.01
            print("[+] mu:",mu)
            result = Coppersmith(f,r,1,beta,mu)
            if result != []:
                return gcd(result[0][0] + (u * inverse_mod(e,N)),N)
            

gamma = 0.06                     # e = N ^ gamma
beta = 0.2                      # Q = N ^ beta
r = 4                             # N = P * Q ^ r
k = 20                           # The number of prime
lb = 20                          # bits-length of primes
Bq = 300                         # bits-length of large prime in P
RR = RealField(2 ** 4)
P,Q,e = Generate(k,lb,Bq,r,beta,gamma)
N = P * Q ** r

gamma = RR(Integer(e).nbits() / Integer(N).nbits())
beta = RR(Integer(Q).nbits() / Integer(N).nbits())

print("[+] beta: {}, gamma: {}".format(beta,gamma))
print("[+] Theo. low bounds of gamma: {}, {}".format( RR(0.25 / r), RR(beta - r * beta ** 2)))
Starttime = time.time()
# print(Theorem1(N,e,r))
print(Theorem1(N,e,r,beta))
print("[+] runtime:",time.time() - Starttime)
