
#######################################################
## Source code for implementations of VSS schemes    ##
## Π_P and Π_LA, built using the general framework   ##
## Π, presented at PKC 2025 (ia.cr/2023/1669). This  ## 
## also includes implementations of the well-known   ##
## Shamir secret sharing scheme, as well as the VSS  ## 
## schemes proposed by Atapoor, Baghery, Cozzo, and  ## 
## Pedersen [ABCP23], and Pedersen [Ped91].          ## 
## Author: Karim Baghery                             ## 
## The source code implementations for VSS schemes   ##
## [Ped91] and [ABCP23] are taken from               ## 
## https://github.com/Baghery/VSS-ABCP23             ##  
#######################################################

from hashlib import sha256,sha1
from time import time

#############################
### ====  ABCP23 VSS  ==== ##
#############################

# global parameters
N = 2^255-19   # Curve 25519
ZN = Integers(N)
LAM = Integers(2^128)
R.<x> = PolynomialRing(ZN)

# routines
def shamir_ABCP(n):
    t = n//2 - 1        # assumes n is even!
    f = R.random_element(degree=t)
    feval = [f(x=i) for i in range(1,n+1)]
    return f,feval

def prover_ABCP(f,feval):
    t = f.degree()
    n = len(feval)
    b = R.random_element(degree=t)
    y = [ [LAM.random_element(),LAM.random_element()] for i in range(n)]
    C,C_ = "",""

    for i in range(1,n+1): # parties are 0..n-1
        bi = b(x=i)
        C += sha256((str(bi)+str(y[i-1][0])).encode()).hexdigest()+str(",")
        C_+= sha256((str(feval[i-1])+str(y[i-1][1])).encode()).hexdigest()+str(",")
    C = C[:-1]
    C_= C_[:-1]
    d = Integer(ZN(int(sha256(str(C+C_).encode()).hexdigest(),16)))
    r = b-d*f
    return C,C_,r,y

def verifier_ABCP(i,C,C_,r,yi,yi_,xi):
    ri = r(x=i)
    d = Integer(ZN(int(sha256(str(C+C_).encode()).hexdigest(),16)))
    Ci  = sha256((str(ri+d*xi)+str(yi)).encode()).hexdigest()
    Ci_ = sha256((str(xi)+str(yi_)).encode()).hexdigest()
    return Ci == C.split(',')[i-1] and Ci_ == C_.split(',')[i-1]

# benchmark function
def benchmark_ABCP23(n):
    N = 2^255-19       # Curve 25519
    ZN = Integers(N)
    LAM = Integers(2^128)
    R.<x> = PolynomialRing(ZN)

    Ts = time()
    f,feval = shamir_ABCP(n)
    Ts = time() - Ts

    Tp = time()
    C,C_,r,y = prover_ABCP(f,feval)
    Tp = time() - Tp

    i = randint(1,n-1) # pick a random verifier

    Tv = time()
    b = verifier_ABCP(i,C,C_,r,y[i-1][0],y[i-1][1],feval[i-1])
    Tv = time() - Tv

    if b == False: print(b)
    return Ts,Tp,Tv


#############################
### ==== Pedersen VSS ==== ##
#############################

# global parameters

p = 2^255-19         # Curve 25519
q = 2^252 + 27742317777372353535851937790883648493 #point group size
Zq = Integers(q)     # group of multiplication map
Am = 486662          # Montgomery A-coefficient
Ar = int((Am+2)/4)   # reduced Montgomery coefficent
E = EllipticCurve(GF(p),[0,Am,0,1,0])
RP.<x> = PolynomialRing(Zq)

G = E.random_point() # generator 1
while G.order() != q:
    G = E.random_point()
s = Zq.random_element()
H = Integer(s)*G     # generator 2

s2 = Zq.random_element()
K = Integer(s2)*G     # generator 3


# Montgomery subroutines

def xADD(P,Q,R): # points are of the form [X,Z]
    [XP,ZP] = [P[0],P[1]];
    [XQ,ZQ] = [Q[0],Q[1]];
    [XR,ZR] = [R[0],R[1]];

    V0 = XP + ZP
    V1 = XQ - ZQ
    V1 = V1 * V0
    V0 = XP - ZP
    V2 = XQ + ZQ
    V2 = V2 * V0
    V3 = V1 + V2
    V3 = V3^2
    V4 = V1 - V2
    V4 = V4^2
    Xp = ZR * V3
    Zp = XR * V4
    
    return [Xp,Zp]

def xDBL(P): # points are of the form [X,Z]
    [XP,ZP] = [P[0],P[1]]
    
    V1 = XP + ZP
    V1 = V1^2
    V2 = XP - ZP
    V2 = V2^2
    X2 = V1 * V2
    V1 = V1 - V2
    V3 = Ar * V1
    V3 = V3 + V2
    Z2 = V1 * V3
    
    return [X2,Z2]

def Montgomery_ladder(k,P): # points are of the form [X,Z]
    x0,x1 = P,xDBL(P)
    k = k.binary()
    l = len(k)
    for i in range(1,l):
        if k[i]=='0':
            x0,x1 = xDBL(x0),xADD(x0,x1,P)
        if k[i]=='1':
            x0,x1 = xADD(x0,x1,P),xDBL(x1)
    return x0,x1

def recover_y(P,Q,R):
    [XP,YP] = [P[0],P[1]] # P is an actual elliptic curve point in the form (X:Y:Z)
    [XQ,ZQ] = [Q[0],Q[1]]
    [XR,ZR] = [R[0],R[1]]
        
    V1 = XP * ZQ
    V2 = XQ + V1
    V3 = XQ - V1
    V3 = V3^2
    V3 = V3 * XR
    V1 = 2*Am*ZQ
    
    V2 = V2 + V1
    V4 = XP * XQ
    V4 = V4 + ZQ
    V2 = V2 * V4
    V1 = V1 * ZQ
    V2 = V2 - V1
    V2 = V2 * ZR
    
    Y  = V2 - V3
    V1 =  2 * YP
    V1 = V1 * ZQ
    V1 = V1 * ZR
    X  = V1 * XQ
    Z  = V1 * ZQ
    
    return E(X,Y,Z)

def fast_multiply(k,P): # use montgomery ladder and y-recovery
    PM = [P[0],P[2]] # X-Z coordinates
    x0,x1 = Montgomery_ladder(Integer(k),PM)
    return E(recover_y(P,x0,x1))


# routines
def shamir_P(n):
    t = n//2-1 # assumes n is even!
    f = RP.random_element(degree=t)
    feval = [f(x=i) for i in range(1,n+1)]
    return f,feval

def prover_P(f,feval):
    t = f.degree()
    C = [0 for _ in range(t+1)]
    g = RP.random_element(degree=t)
    fcoeff = f.coefficients(sparse=False)
    gcoeff = g.coefficients(sparse=False)
    geval = [g(x=i) for i in range(1,n+1)]
    for i in range(t+1):
        C[i] = fast_multiply(Integer(fcoeff[i]),G) + fast_multiply(Integer(gcoeff[i]),H)
    return g,C

def verifier_P(a,b,C,i):
    S = C[0]
    for j in range(1,len(C)):
        S += fast_multiply(Integer(Zq(i)^j),C[j])
    return S == fast_multiply(a,G) + fast_multiply(b,H)


# benchmark function
def benchmark_Pedersen(n):
    Ts = time()
    f,feval = shamir_P(n)
    Ts = time() - Ts

    Tp = time()
    g,C = prover_P(f,feval)
    Tp = time()-Tp

    i = randint(1,n) # sample verifier
    a,b = Integer(f(x=i)),Integer(g(x=i))

    Tv = time()
    boo = verifier_P(a,b,C,i)
    Tv = time()-Tv

    if boo == False: print(boo)
    return Ts,Tp,Tv


#############################
##### ==== Pi_P VSS ==== ####
#############################

# routines
def shamir_Pi_P(n):
    t = n//2-1 # assumes n is even!
    f = RP.random_element(degree=t)
    feval = [f(x=i) for i in range(1,n+1)]
    return f,feval

def prover_Pi_P(f,feval):
    t = f.degree()
    C = [0 for _ in range(n)]
    g = RP.random_element(degree=t)
    geval = [g(x=i) for i in range(1,n+1)]
    gamma_rand = [Zq.random_element() for _ in range(1, n+1)]   # randomizers for the commitments 
    for i in range(n):
        C[i] = fast_multiply(Integer(feval[i]),G) + fast_multiply(Integer(geval[i]),H) + fast_multiply(Integer(gamma_rand[i]),K)   
    d = Integer(Zq(int(sha256(str(C).encode()).hexdigest(),16)))
    z = g+d*f       
    return g,C, gamma_rand, z

def verifier_Pi_P(i,Ci,C,z,fi,gammai):
    d = Integer(Zq(int(sha256(str(C).encode()).hexdigest(),16)))
    zi = Integer(z(x=i))    
    gi = Integer(zi - d*fi) 
    return Ci == fast_multiply(fi,G) + fast_multiply(gi,H) + fast_multiply(gammai,K)


# benchmark function
def benchmark_Pi_P(n):
    Ts = time()
    f,feval = shamir_Pi_P(n)
    Ts = time() - Ts

    Tp = time()
    g,C, gamma_rand, z = prover_Pi_P(f,feval)
    Tp = time()-Tp

    i = randint(1,n) # sample verifier
        
    Ci = C[i-1]    # commit to shares of randomizer 
    fi = f(x=i)
    gammai = Integer(gamma_rand[i-1])  # randomizer of party i in commitment Ci  
    Tv = time()
    boo = verifier_Pi_P(i,Ci,C,z,fi,gammai)
    Tv = time()-Tv

    if boo == False: print(boo)    
    return Ts,Tp,Tv


################################################
## Pi_LA VSS (with high-entropy secret, i.e., ## 
##            without the randomizer \gamma)  ##
################################################

def shamir_LA(n):
    t = n//2 - 1 # assumes n is even!
    f = RP.random_element(degree=t)
    feval = [f(x=i) for i in range(1,n+1)]
    return f,feval

def prover_LA(f,feval):
    t = f.degree()
    n = len(feval)
    b = RP.random_element(degree=t)
    C = ""

    for i in range(1,n+1): # parties are 0..n-1
        bi = b(x=i)
        C+= sha256((str(feval[i-1]) + str(bi)).encode()).hexdigest()+str(",")
    C= C[:-1]
    d = Integer(Zq(int(sha256(str(C).encode()).hexdigest(),16)))
    r = b-d*f    
    return C,r

def verifier_LA(i,C,r,xi):
    ri = r(x=i)
    d = Integer(Zq(int(sha256(str(C).encode()).hexdigest(),16)))
    Ci = sha256((str(xi)+str(ri+d*xi)).encode()).hexdigest()
    return Ci == C.split(',')[i-1]

# benchmark function
def benchmark_LA(n):
    N = 2^256-189
    ZN = Integers(N)
    LAM = Integers(2^128)
    R.<x> = PolynomialRing(ZN)

    Ts = time()
    f,feval = shamir_LA(n)
    Ts = time() - Ts

    Tp = time()
    C,r = prover_LA(f,feval)
    Tp = time() - Tp

    i = randint(1,n-1) # pick a random verifier

    Tv = time()
    b = verifier_LA(i,C,r,feval[i-1])
    Tv = time() - Tv

    if b == False: print(b)
    return Ts,Tp,Tv


# running the benchmark
N=[16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
# N=[32, 128, 512, 2048, 8192, 16384]   # Reported in the paper 

for i in range(len(N)):
  n = N[i]  
  iteration = 1
  if n < 513: iteration = 10

  print("=========================================================================================")
  print("Benchmarking Shamir and VSS Schemes Pedersen, Pi_P, ABCP, Pi_LA for (n, t) =", (n, n/2-1))
  print("=========================================================================================")

  Ts1f = 0 
  Tp1f = 0
  Tv1f = 0
  
  Ts2f = 0 
  Tp2f = 0
  Tv2f = 0
  
  Ts3f = 0 
  Tp3f = 0
  Tv3f = 0  

  Ts4f = 0 
  Tp4f = 0
  Tv4f = 0
  
  # Ts,Tp,Tv = benchmark_Pedersen(n)
  # print("Pedersen scheme (public verification) -- proof time:       ",Ts+Tp)
  # print("Pedersen scheme (public verification) -- verification time:",Tv)
  for i in range(iteration):
    Ts1,Tp1,Tv1 = benchmark_Pedersen(n)
    Ts2,Tp2,Tv2 = benchmark_ABCP23(n)
    Ts3,Tp3,Tv3 = benchmark_LA(n)
    Ts4,Tp4,Tv4 = benchmark_Pi_P(n) 
    Ts1f += Ts1 
    Tp1f += Tp1
    Tv1f += Tv1
    
    Ts2f += Ts2 
    Tp2f += Tp2
    Tv2f += Tv2
    
    Ts3f += Ts3 
    Tp3f += Tp3
    Tv3f += Tv3
    
    Ts4f += Ts4 
    Tp4f += Tp4
    Tv4f += Tv4
    

  print("       ======================== Shamir Secret Sharing ===========================")
  print("       Shamir    -- sharing time:       ",(Ts2f)/iteration)
  print("       ==========================================================================")

  print("       ======================== VSS Schemes - Sharing Phase =====================")
  print("       Pedersen  -- sharing time:       ",(Ts1f+Tp1f)/iteration)
  print("       Pi_P      -- sharing time:       ",(Ts4f+Tp4f)/iteration)
  print("       ABCP23    -- sharing time:       ",(Ts2f+Tp2f)/iteration)
  print("       Pi_La     -- sharing time:       ",(Ts3f+Tp3f)/iteration)

  print("       ======================== VSS Schemes - Verification ======================")
  print("       Pedersen  -- verification time:  ",Tv1f/iteration)
  print("       Pi_P      -- verification time:  ",Tv4f/iteration)
  print("       ABCP23    -- verification time:  ",Tv2f/iteration)
  print("       Pi_La     -- verification time:  ",Tv3f/iteration)
  print("       ==========================================================================")