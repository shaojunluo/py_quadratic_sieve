from math import log

# ==================Helper functions================#

# Helper function to fast evaluate a^b mod m
def power_mod(a,b,m):
    n = 1
    while( b!= 0):
        if( b % 2 == 1 ):
            n = n*a % m
        b = b // 2
        a = a*a % m
    return n

# find a non-residue of p
def find_non_residue(p):
    for b in range(p):
        if(jacobi_symbol(b,p) == -1):
            return b

# helper function for transposing a matrix
def transpose(A):
    A_T = []
    for i in range(len(A[0])):
        row_T = []
        for row in A:
            row_T.append(row[i])
        A_T.append(row_T)
    return A_T

# helper function for reading files
def read_file(file_name):
    with open(file_name,'r') as f:
        digits = f.readlines()[0].strip().replace('[','').replace(']','')
        digits = digits.split(',')
    M,row = [],[] # construct matrix
    for i,x in enumerate(digits):
        row.append(int(x))
        if (i+1) % 100 == 0:
            M.append(row)
            row = []
    return M

# vector power to construct number
def vec_pow(a, s):
    P = 1
    for i in range(len(a)):
        P *= a[i]**s[i]
    return P

# helper function for integer sqrt use newton's method
def isqrt(n):
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

# helper function to compute gcd
def gcd(x, y): 
    while(y): 
        x, y = y, x % y 
    return x 

# ===================Problem 1 =====================#

# P1(a) Factorize given a list of prims
def factorize_with_primes(n,ps):
    indices = []
    for p in ps:
        idx = 0
        while n % p ==0:
            n = n // p
            idx += 1
        indices.append(idx)
    return indices, n

# P1(b) pseudoprime_test with Miller-Rabin
def pseudoprime_test(n, trials):
    # some null cases
    if n ==1:
        return False
    elif n == 2:
        return True
    elif n % 2 == 0:
        return False
    
    # Factorize n-1 as 2^s*d
    s,d = factorize_with_primes(n-1,[2])
    s = s[0]
    
    # search for witness
    for a in trials:
        is_witness = True
        # condition 1
        if power_mod(a, d, n) == 1:
            is_witness = False
            continue
        # condition 2
        for i in range(s):
            if power_mod(a, 2**i*d, n) == n-1:
                is_witness = False
                break
        
        # if both condition not satisfied, a is witness
        if is_witness:
            print(f'Target: {n} Pseudo Test Result: witness_found: {a}, Not prime')
            return False
        
    # else cannot assert it is composite
    print(f'n = {n}, Pseudo Test Result: NO witness_found, Probaly Prime')
    return True

# ======================= P2 ========================#

# P2(a) Code for jacobi symbol from HW 10
def jacobi_symbol(n, m):
    assert m%2==1
    ans = 1
    while n > 1:
        # t is the number of 2 pulled from n
        t = 0
        while n % 2 == 0:
            n = n // 2
            t += 1
        # get the result Jacobi_symbol(2, m) pow t times
        if ((m * m - 1) / 8) % 2 == 1 and t % 2 == 1:
            ans = -1*ans
        if n == 1:
            break
        #swap and reduce n,m
        if n >= m:
            n = n % m
        else:
            if (((n-1)*(m-1))//4) % 2 == 1:
                ans = -1 *ans
            # swap m,n
            temp = m
            m = n
            n = temp
    if n == 0:
        return 0
    return ans

# P2 (b) Extended Euclidean algorithm (a la Knuth; see book)
def extended_euclid(m,n):
    olda, a = 1, 0
    oldb, b = 0, 1
    while n > 0:
        q = m // n
        a, olda = olda - q*a, a
        b, oldb = oldb - q*b, b
        m, n = n, m % n
    return m, olda, oldb

# P2 (d） Tonelli Algorithm
def tonelli(a,p,verbose= 0):
    if verbose:
        print(f"Tonelli Algorithm for square root of {a} modulo {p}",end = '...')
    if(jacobi_symbol(a,p) != 1 ):
        if verbose:
            print(f"There doesin exit a squared number with {a} mod {p} , quit")
        return -1, []
    s,d = factorize_with_primes(p-1,[2])
    s = s[0]
  
    # calculate inverse modulos
    b = find_non_residue(p)
    _, binv, _ = extended_euclid(b,p)
    binv = binv % p
    
    # Search for i sp that (a/b^i)^t = 1 mod p")
    i_s = [0,2]
    # index the i_s so i[j] = i_j by sticking a dummy number in i[0].
    for k in range(2,s+1):
        check = power_mod(a*power_mod(binv,i_s[k-1],p),d*(2**(s-k)),p)
        if(check == 1):
            i_s.append(i_s[k-1])
        else:
            i_s.append(i_s[k-1]+2**(k-1))
            
    root = power_mod(b,i_s[-1]/2,p) * power_mod(a*power_mod(binv,i_s[-1],p),(d+1)/2,p) % p
    if verbose:
        print(f"Sqrt Root found: {root}")
    return root, i_s

#====================== Problem 3================================#

# P3(a) Faset Gaussian elimination under binary matrix (modulus 2). Implementation of:
# https://www.cs.umd.edu/~gasarch/TOPICS/factoring/fastgauss.pdf
def gaussian_elimination_mod2(M):

    contains_pivot = [False]*len(M)
    pivot_found = False
    pivot_loc = {}
    
    for j in range(len(M[0])):
        pivot_found = False
        #Look for pivot
        for i in range(len(M)):
            #Pivot Found at row i and column j
            if(M[i][j] == 1):
                contains_pivot[i] = True
                pivot_loc[j]=i
                pivot_found = True
                break
          
        if pivot_found:
            for k in range(len(M[0])):
                if (k == j):
                    continue
                if (M[i][k] == 1):
                    for l in range(len(M)):
                        M[l][k] = (M[l][j] + M[l][k])%2
                        
    return M, contains_pivot, pivot_loc

# P3(b) find the left null space of M mod 2
def left_null_space_mod2(M): 
    # solve with gaussian elimination
    M, contains_pivot, pivot_loc = gaussian_elimination_mod2(M)
    result = []
    # Find left null space
    for i in range(len(M)):
        #Find dependent rows
        null_space = []
        if not contains_pivot[i]:
            row = M[i]
            null_space = [i]
            for j in range(len(row)):
                if (M[i][j]==1):
                    null_space.append(pivot_loc[j])
            result.append(null_space)
    # Encoding the indicies
    out_vecs = []
    for space in result:
        vec = [0]*len(M)
        for idx in space:
            vec[idx] = 1
        out_vecs.append(vec)
    return out_vecs
            
#====================== Problem 4 ================================#

# P4(a,b) List all primes p that smaller than B and Jacobi symbol(n,p) = 1. 
# If n is 1 then list all primes.
def list_primes(B, n = 1):
    if B < 2:
        return []
    if B == 2:
        return [2]
    
    # use compact sieve of eratosthenes 
    sieve = [True] * (B//2)
    i = 3
    while i*i <= B:
        # mask seive
        if sieve[i//2]:
            sieve[i*i//2::i] = [False]*((B-i*i-1)//(2*i)+1)
        i += 2
    
    p_list = []
    for p in range(1, B//2):
        if sieve[p] and (jacobi_symbol(n,2*p+1)==1):
            p_list.append(2*p+1)
      
    return [2] + p_list

#====================== Problem 5 ================================#

# P5(a) generate factor base
def calculate_h(n, M):
    start = int(n**0.5)+1
    hs = []
    rs = list(range(start, start + M +1))
    for r in rs:
        hs.append(log(r*r -n))
    return hs, rs

# P5(b) generate factor base
def solve_factor_base(n, B):
    tps = {}
    ps = list_primes(B, n = n)
    for p in ps[1:]:
        root, _ = tonelli(n,p)
        tps[p] = root
    return tps, ps

# P5(de) run sieve. We replace r with t if t is smaller
def run_sieve(n,M,B,tol = 1e-3):
    #  null cases
    if n < 2:
        return []
    if n == 2:
        return [2]
    
    # value of h and t
    hs,rs = calculate_h(n,M)
    ts,ps = solve_factor_base(n,B)
    factor_idx,x,y  = [],[],[]
    
    # Sieving
    Y = [(r*r - n) for r in rs]
    for j,p in enumerate(ps):
        # Perform sieve to reduce hs value
        for i in range(len(Y)):
            if Y[i] % p == 0:
                hs[i] -= log(p)
                Y[i] = Y[i]//p
                
    # search candidate
    for i in range(len(hs)):
        # find factors in rs
        if abs(hs[i]) < tol:
            s,d = factorize_with_primes(rs[i]*rs[i]-n,ps)
            factor_idx.append(s)
            # replace r as t values to reduce size
            x_val = vec_pow(list(ts.values()),s)
            if (rs[i]*rs[i]-n) % 2 == 0:
                x_val *= 2
            x.append(min(rs[i],x_val))
            y.append(min(abs(rs[i]*rs[i]-n),abs(x_val*x_val-n)))
    return factor_idx,x,y

# =====================Problem 6=======================#

def quadratic_sieve(n,M,B,trials):
    pseudoprime_test(n, trials)
    # RUN sieve to get matrix
    print(f'Running sieve with M = {M}, B = {B}')
    A,x,y = run_sieve(n,M,B,tol=0.01)
    # find the solution of indicies
    print('Solve left null space using Gaussian on GF(2)...')
    coefs = left_null_space_mod2(A)
    x = [vec_pow(x, c) for c in coefs]
    y = [vec_pow(y, c) for c in coefs]
    print(f'{len(x)} candidates found. Searching...')
    for k in range(len(x)):
        factor = gcd(isqrt(y[k])-x[k],n)
        if (factor != n) & (factor != 1):
            print(f'Factor found: {factor}')
            print(f'x = {x[k]} \ny = {isqrt(y[k])}')
            return factor
    print(f'NO factors found, try increasing B or M')
    