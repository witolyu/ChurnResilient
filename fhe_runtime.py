import math

from scipy.stats import binom

n = 10000
t = 3/5
w = 10000
d = 10000
lam = 40
deg = 2**14

m_range = [50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000]
m_polynomial_times= [0.0093226,0.0133125,0.016475,0.0216151,0.0286177,0.0265638,0.0300709,0.0362433,0.0428097,0.0487032,
0.0570238,0.0630198,0.0664729,0.0598916,0.0624344,0.066745,0.070887,0.0775334,0.0988227,0.0976743]

for i,m in enumerate(m_range):
    # m is the size of decryption committee.
    # where mc is the size of computation committee.

    # We assert there is enough parties so that every computation committee has their own decrypotion committee.
    # number of committee

    c = math.ceil(w / deg)
    assert m*c < n


    # this is the "negligible" probability that each decryption committee fail to contains the number of corruption.
    op = 1 - 2 ** (-lam - 1)/c

    # Initial corrupt parties, at least larger than expected value.
    cor = math.ceil(m * t)

    while binom.cdf(cor,m,t) < op:
        cor += 1
        if cor == m-1:
            break
    print("-----------")
    if cor == m-1:
        print("Committee size",m,"not large enough")
    else:
        # note that we should set cor+1 as the threshold. (Use cor-degree polynomial), this allows m-cor -1 failures.
        # However, we would like to test different fraction of the maximum failure parties, as we assume pessimisitic
        # that all of them will fail. This gives us the amoumt of reconstruction that needs to be done.

        for f in [int((m-cor)/8),int((m-cor)/4),int((m-cor)/8*3),int((m-cor)/2),int((m-cor)/4*3)]:
            # this is the time needed for reconstruction.
            tau = f * m_polynomial_times[i]

            # add the threshold decryption time, is the computational time independent of the committee size?
            decryption_time = 0.1
            tau_w_dec = tau + d* 0.1

            # this computes the tolerable end-to-end failure rate per party.
            p = 0
            while binom.cdf(f,m,p) > 1-op:
                p += 0.01

            alpha = 1-(1-p)**(1/tau_w_dec)
            # Printing variables with 6 slots allocated for each
            print(f"Committee size: {m:<6}, Corrupted parties: {cor:<6}, Allowed failed parties: {f:<6}, "
                  f"recovery time(Interpolation): {tau:<6}, decrpytion time {tau_w_dec-tau:<6}," 
                  f"Decryption committee runtime: {tau_w_dec:<6}, Tolerable failure rate per second:", alpha)






