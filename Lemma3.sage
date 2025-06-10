from tqdm import tqdm
from re import findall
from subprocess import check_output
def flatter(M):
    # compile GitHub - keeganryan/flatter: Fast lattice reduction and put it in $PATH
    z = "[[" + "]\n[".join(" ".join(map(str, row)) for row in M) + "]]"
    ret = check_output(["flatter"], input=z.encode())
    return matrix(M.nrows(), M.ncols(), map(int, findall(b"-?\\d+", ret)))



def find_roots(B,method = 'variety',bound = None,monomials = None,f = None):
    PR.<x> = PolynomialRing(ZZ)
    roots = []
    if method == 'variety':
        H = Sequence([], f.parent().change_ring(QQ))
        for h in filter(None, B*monomials):
            H.append(h)
            I = H.ideal()
            if I.dimension() == -1:
                H.pop()
            elif I.dimension() == 0:
                for root in I.variety(ring=ZZ):
                    root = tuple(root[var] for var in [x])
                    roots.append(root)
    elif method == 'sage':
        h = sum([ZZ(B[0,-(i+1)]) * x**i for i in range(B.ncols())])
        roots = h.roots()
    return roots

def Coppersmith(f,u,v,beta,epsilon):
    m = ceil((beta * (2 * u + v - u * v * beta))/(epsilon)) - 1
    t = ceil(u * beta * m)
    # print(f"[+] m = {m}, t = {t}")
    R = f.base_ring()
    N = R.cardinality()
    f /= f.coefficients().pop(0)
    f = f.change_ring(ZZ)
    bounds = [int(N ^ (u * v * beta ^ 2 - epsilon))]

    g = lambda k: f ^ k * N ^ max(ceil((v * (t - k)) / u ),0)
    ShiftPolys = Sequence([],f.parent())
    for k in range(m + 1):
        Poly = g(k)
        ShiftPolys.append(Poly)
    B, monomials = ShiftPolys.coefficients_monomials()
    monomials = vector(monomials)
    nn = len(ShiftPolys)
    print("[+] dim:",B.dimensions())

    factors = [monomial(*bounds) for monomial in monomials]
    for i,factor in enumerate(factors):
        B.rescale_col(i,factor)

    # B = B.dense_matrix().LLL(delta = 0.75)
    
    B = flatter(B)
    # print('[+] LLL done')
    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1/factor)

    # roots = find_roots(B,method = 'variety', monomials = monomials,f = f)
    roots = find_roots(B,method = 'sage')
    return roots

def iterationAlgorithm(f,r,gamma,mu):
    print(f'[+] Start iterationAlgorithm.')
    beta0 = 1 / r
    Nextbeta = lambda beta: max(sqrt((beta - (gamma - mu))/ r),gamma)
    betas = [beta0]

    while True:
        betan = Nextbeta(beta0)
        beta0 = betan
        if beta0 == gamma:
            break
        betas.append(beta0)
    print(f'[+] Total {len(betas)} steps. The max of m: {ceil((betas[0] * (2 * r + 1 - r * 1 * betas[0] ))/(mu)) - 1}')
    roots = []
    for beta in tqdm(betas):
        res = Coppersmith(f,r,1,beta,mu)
        res = list(set([int(root[0]) for root in res]))
        for root in res:
            roots.append(root)
    return roots
