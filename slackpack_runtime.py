import math

from scipy.stats import binom

n = 10000
t = 3/5
w = 10000
d = 10000
lam = 40

m_range = [50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000]
m_polynomial_times= [0.0093226,0.0133125,0.016475,0.0216151,0.0286177,0.0265638,0.0300709,0.0362433,0.0428097,0.0487032,
0.0570238,0.0630198,0.0664729,0.0598916,0.0624344,0.066745,0.070887,0.0775334,0.0988227,0.0976743]

for i,m in enumerate(m_range):
    # number of committee
    c = math.floor(n / m)

    # compute number of corrupt parties per committees
    cor = math.ceil(m*t)

    op = 1 - 2 ** (-lam - 1)/c

    while binom.cdf(cor,m,t) < op:
        cor += 1
        if cor == m-1:
            break
    print("-----------")
    if cor == m-1:
        print("Committee size",m,"not large enough")
    else:
        for f in [int((m-cor)/8),int((m-cor)/4),int((m-cor)/8*3),int((m-cor)/2),int((m-cor)/4*3)]:
            gamma = int((m-cor - f +1)/2)
            l = math.ceil(w/c/gamma)
            tau = d * l * m_polynomial_times[i]

            p = 0
            while binom.cdf(f,m,p) > 1-op:
                p += 0.01

            alpha = 1-(1-p)**(1/tau)
            # Printing variables with 6 slots allocated for each
            print(f"Committee size: {m:<6}, Corrupted parties: {cor:<6}, Slack parties: {f:<6}, Packing parameters: {gamma:<6}, pack/committee {gamma/m:<6}," 
                  f" Batch size: {l:<6}, Eval/Inter Runtime: {int(tau):<6}, Tolerable failure rate per second:", alpha)






