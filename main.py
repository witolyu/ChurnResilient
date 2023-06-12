# This is a sample Python script.
import math
from scipy.stats import binom
import matplotlib.pyplot as plt


# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def CommiteeSize(fail_rate_per_ms, runTime, t, lam):
    # Assuming selection with replacement.
    # using union bound to compute the fail rate? Is it too loose?
    # Compute the size needed to ensure one honest survival with 2^-lam probability.

    if fail_rate_per_ms * runTime >= 1:
        return -1

    bad = t + (1 - t) * (runTime * fail_rate_per_ms)

    fail_committee_rate = bad
    commiteeSize = 1
    while fail_committee_rate > 2 ** (-40):
        commiteeSize += 1
        fail_committee_rate *= bad

    return commiteeSize


def CommiteeSizeFullPartition(n, t, fail_rate_per_ms, runTime, lam):
    # Union bound

    # Partition
    return -1


def FullContriRunTime(n, t, scalarTime, B):
    # compute the runtime to do a vandermonde matrix multiplication of size ((1-t)u * u)

    # FIXME: we probably also need to add ciphertext multiplication

    return (1 - t) * n * n * scalarTime * B


def FullContriTriples(n, t, m, l, B):
    # compute the number of triples generate

    # FIXME: adjust the triple number apart from random values number.

    return (1 - t) * n * l * m * B


def binary_search_largest_smaller(left, right, target, u, p):
    largest_smaller = None
    while left <= right:
        mid = (left + right) // 2
        if binom.cdf(mid, u, p) < target:
            largest_smaller = mid
            left = mid + 1
        else:
            right = mid - 1
    return largest_smaller


def PartContriRunTime(u, t, scalarTime, B, lam):
    # compute the honest contributor lower bound.
    h = binary_search_largest_smaller(0, math.floor(u * (1 - t)), 2 ** (-lam), u, (1 - t))

    # compute the runtime to do a vandermonde matrix multiplication of size ((1-t)u * u)

    # FIXME: we probably also need to add ciphertext multiplication and Zero Knowledge Proof.

    return h * u * scalarTime * B


def PartContriTriples(u, t, m, l, B, lam):
    # compute the honest contributor lower bound.
    h = binary_search_largest_smaller(0, math.floor(u * (1 - t)), 2 ** (-lam), u, (1 - t))
    # FIXME: adjust the triple number apart from random values number.

    return h * l * m * B


def Test1():
    lam = 40

    t = 2.0 / 3  # corruption threshold
    n = 10000  # network size
    scalarTime = 11  # ms to compute a ciphertext plaintext multiplication.
    l = 2 ** 15  # packing parameter

    # Case 1, Full network contribution, full network participation, 1 batch.
    B = 1
    runTime = FullContriRunTime(n, t, scalarTime, 1)
    fail_rate_per_ms = 0.1 / runTime
    c = CommiteeSize(fail_rate_per_ms, runTime, t, lam)
    m = math.floor(n / c)
    T = FullContriTriples(n, t, m, l, B)
    print(("Case 1: runTime = {}ms, fail_rate_per_ms={}, c = {},  m = {}, T = {}".format(runTime, fail_rate_per_ms, c,
                                                                                         m, T)))
    print("1.0", runTime, fail_rate_per_ms, c, m, T)

    # Case 2, Full network contribution, full network participation, 2 batch.
    B = 2
    runTime = FullContriRunTime(n, t, scalarTime, 2)
    # fail_rate_per_ms = 0.1 / runTime
    c = CommiteeSize(fail_rate_per_ms, runTime, t, lam)
    m = math.floor(n / c)
    T = FullContriTriples(n, t, m, l, B)
    # print(("Case 2: runTime = {}ms, fail_rate_per_ms={}, c = {},  m = {}, T = {}".format(runTime, fail_rate_per_ms, c,
    #                                                                                      m, T)))

    # Case 3, Partial network contribution, full network participation, 1 batch.

    for ratio in range(9, 0, -1):
        B = 1
        u = n * ratio / 10
        runTime = PartContriRunTime(u, t, scalarTime, 1, lam)
        c = CommiteeSize(fail_rate_per_ms, runTime, t, lam)
        m = math.floor(n / c)
        T = PartContriTriples(u, t, m, l, B, lam)
        # print(("Case 3: contribution_ratio={}, runTime = {}ms, fail_rate_per_ms={}, c = {},  m = {}, T = {}".format(ratio/10, runTime, fail_rate_per_ms, c,
        #                                                                                      m, T)))
        print(ratio / 10, runTime, fail_rate_per_ms, c, m, T)


# FIXME All the functions above may not be necessary, consider deleting them once push to github.

plotGraph = 1

def ComputeNpiFloat(t, c, rho):
    # This computes the float n' needed in order to guarantee c ciphertexts.
    lnp = math.log(rho)
    npi = (- lnp + c - 1 + math.sqrt(lnp ** 2 - 2 * lnp * (c - 1))) / (1 - t)
    assert npi * (1 - t) >= c
    return npi


def ComputeNpi(t, c, rho):
    # This computes the n' needed in order to guarantee c ciphertexts.
    return math.ceil(ComputeNpiFloat(t, c, rho))


def Compute_delta1(npi, t, c):
    return math.e ** (-(1 - (c - 1) / (1 - t) / npi) ** 2 * (1 - t) * npi / 2)


def ComputeAlpha(n, t, m, npi, c, rho):
    assert n * t > m
    return 1 - (1 - (rho ** (1 / m) - n * t / (n - m)) * n / (n - n * t + m)) ** (1 / npi / c)


def Compute_delta2(n, t, m, npi, c, alpha):
    p = 1 - (1 - alpha) ** (npi * c)
    return (n * t / (n - m) + (1 - (n * t - m) / n) * p) ** m


def ComputeRuntime(npi, c, beta):
    # This computes the vandermonde runtime,
    return npi * c * beta


def ComputeRuntimeAndAlpha(n, t, circuit_size, d, k, beta, kappa):
    c = math.ceil(circuit_size / d / k / 2)
    assert c < n
    rho = 2 ** (-kappa) / 2 / k
    npi = ComputeNpi(t, c, rho)
    m = math.floor(n / k)
    alpha = ComputeAlpha(n, t, m, npi, c, rho)
    runtime = ComputeRuntime(npi, c, beta)
    return alpha, runtime


def ComputeRuntimesAndAlphas(n, t, circuit_size, d, beta, kappa):
    res, x, y1, y2 = [],[],[],[]

    for k in range(20,n,5):
        c = math.ceil(circuit_size / d / k / 2)
        assert c < n
        rho = 2 ** (-kappa) / 2 / k
        npi = ComputeNpi(t, c, rho)
        m = math.floor(n / k)
        if n * t <= m:
            continue
        if t**m > rho:
            break
        alpha = ComputeAlpha(n, t, m, npi, c, rho)
        if alpha < 0:
            continue
        runtime = ComputeRuntime(npi, c, beta)
        if plotGraph:
            x.append(k)
            y1.append(alpha)
            y2.append(runtime/1000) # Measured in seconds
        res.append((k, m, npi,alpha, runtime))

    if plotGraph:
        plotGraphForAlphasAndRuntimes(x,y1,y2)
    return res

def plotGraphForAlphasAndRuntimes(x, y1, y2):
    # Create figure and axes objects
    fig, ax1 = plt.subplots()

    # Plot data on the first axis
    ax1.plot(x, y1, 'b-', label='alpha')
    ax1.set_xlabel('Number of committee')
    ax1.set_ylabel('tolerable failure rate per scalar multiplication', color='b')
    ax1.tick_params('y', colors='b')
    ax1.minorticks_on()
    # ax1.set_yscale('log')  # Set logarithmic scale on y-axis
    

    # Create a second axis sharing the same x-axis
    ax2 = ax1.twinx()

    # Plot data on the second axis
    ax2.plot(x, y2, 'r-', label='runtime (s)')
    ax2.set_ylabel('runtime (s)', color='r')
    ax2.tick_params('y', colors='r')
    ax2.minorticks_on()

    # ax1.grid(True, which='both')
    # ax2.grid(True, which='both')
    # ax2.set_ylim(ax1.get_ylim())
    # Add legend and title
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    plt.title('Tradeoff between tolerable per scalar multiplication failure rate and runtime')

    # Display the plot
    plt.show()

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('Churn Resilient')

    n = 10000
    t = 2 / 3
    # k = 100
    # m = math.floor(n / k)
    # c = 50

    circuit_size = 100000000
    d = 10000
    # 9 ms per plaintext ciphertext multiplication.
    beta = 9
    kappa = 40

    # npi = ComputeNpi(t, c, 2 ** (-41))
    #
    # alpha = ComputeAlpha(n, t, m, npi, c, 2 ** (-41))
    #
    # print(npi, alpha)

    ComputeRuntimesAndAlphas(n, t, circuit_size, d, beta, kappa)


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
