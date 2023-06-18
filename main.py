# This is a sample Python script.
import math
from scipy.stats import binom
import matplotlib.pyplot as plt
import os


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
saveGraph = 1
folder_path = "Images/"
if not os.path.exists(folder_path):
    os.makedirs(folder_path)


# \label{col:vdm_extract}
def ComputeNpiFloat(t, c, rho):
    # This computes the float n' needed in order to guarantee c ciphertexts.
    lnp = math.log(rho)
    npi = (- lnp + c - 1 + math.sqrt(lnp ** 2 - 2 * lnp * (c - 1))) / (1 - t)
    assert npi * (1 - t) >= c
    return npi


def ComputeNpi(t, c, rho):
    # This computes the n' needed in order to guarantee c ciphertexts.
    return math.ceil(ComputeNpiFloat(t, c, rho))


# \label{lem:vdm_extract}
def Compute_delta1(npi, t, c):
    return math.e ** (-(1 - (c - 1) / (1 - t) / npi) ** 2 * (1 - t) * npi / 2)


# \label{lem:vdm_runtime}
def ComputeRuntimeSecond(l, npi, c, beta1, beta2, beta3, beta4, beta5, beta6):
    # This computes the vandermonde runtime (ms) for l matices of size npi * c,
    return (l * npi * (beta1 + beta2) +
            l * npi * c * (beta3 + beta4) +
            l * c * (2 * beta5 + 3 * beta6)) / 1000


# \label{col:committee_fail}
def ComputeAlpha(n, t, m, v, rho):
    assert n * t > m
    return 1 - (1 - (rho ** (1 / m) - n * t / (n - m)) * n / (n - n * t + m)) ** (1 / v)


# \label{lem:committee_fail}
def Compute_delta2(n, t, m, v, alpha):
    p = 1 - (1 - alpha) ** v
    return (n * t / (n - m) + (1 - (n * t - m) / n) * p) ** m


# \label(lem:vdm_runtime{}
def ComputeRuntimeAndAlpha(n, t, circuit_size, k, l, d,
                           beta1, beta2, beta3, beta4, beta5, beta6, kappa):
    c = math.ceil(circuit_size / k / l / d * 2)
    assert c < n
    rho = 2 ** (-kappa) / 2 / k
    npi = ComputeNpi(t, c, rho)
    m = math.floor(n / k)
    runtime = ComputeRuntimeSecond(l, npi, c, beta1, beta2, beta3, beta4, beta5, beta6)
    alpha = ComputeAlpha(n, t, m, runtime, rho)
    return alpha, runtime

# Does not use Vandermonde
def ComputeRuntimeAndAlphaBase(n, t, circuit_size, k, d,
                           beta1, beta2, beta4, beta5, beta6, kappa):
    c = 1
    l = math.ceil(circuit_size / k / d * 2)
    assert l < n
    rho = 2 ** (-kappa) / 2 / k
    npi = ComputeNpi(t, c, rho)
    m = math.floor(n / k)
    runtime = ComputeRuntimeSecond(l, npi, c, beta1, beta2, 0, beta4, beta5, beta6)
    alpha = ComputeAlpha(n, t, m, runtime, rho)
    return alpha, runtime

# Fixed iteration l and test different committee size.
def ComputeRuntimesAndAlphasForK(n, t, circuit_size, l, d,
                             beta1, beta2, beta3, beta4, beta5, beta6, kappa):
    res, x, y1, y2 = [], [], [], []

    for k in range(20, n, 5):
        c = math.ceil(circuit_size / k / l / d * 2)
        assert c < n
        rho = 2 ** (-kappa) / 2 / k
        npi = ComputeNpi(t, c, rho)
        m = math.floor(n / k)
        if n * t <= m:
            continue
        if t ** m > rho:
            break

        runtime = ComputeRuntimeSecond(l, npi, c, beta1, beta2, beta3, beta4, beta5, beta6)
        alpha = ComputeAlpha(n, t, m, runtime, rho)
        if alpha < 0:
            continue

        if plotGraph:
            x.append(k)
            y1.append(alpha)
            y2.append(runtime)  # Measured in seconds
        res.append((k, m, npi, alpha, runtime))

    if plotGraph:
        if saveGraph:
            filename = folder_path + f"RuntimeAlphaVSK_l={l}.png"
            plotGraphForAlphasAndRuntimes(x, y1, y2, filename)
        else:
            plotGraphForAlphasAndRuntimes(x, y1, y2)
    return res

def plotGraphForAlphasAndRuntimes(x, y1, y2, filename = ''):
    # Create figure and axes objects
    fig, ax1 = plt.subplots()

    # Plot data on the first axis
    ax1.plot(x, y1, 'b-', label='alpha')
    ax1.set_xlabel('Number of committee')
    ax1.set_ylabel('tolerable failure rate per second', color='b')
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
    plt.title(f'Tradeoff between tolerable failure rate per second and runtime')

    # Display the plot
    if filename == '':
        plt.show()
    else:
        plt.savefig(filename)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('Churn Resilient')

    n = 10000
    t = 2 / 3
    # k = 100
    # m = math.floor(n / k)
    # c = 50

    circuit_size = 10000000
    d = 2**13

    # The following specifies the runtime (ms) for each type of operation

    # β1 Receive ciphertext
    beta1 = 0

    # β2 Verify zero-knowledge proof
    beta2 = 0

    # β3 Scalar multiplication
    beta3 = 9

    # β4 Ciphertext addition
    beta4 = 0

    # β5 Ciphertext multiplication
    beta5 = 34

    # β6 Unpacking ciphertext
    beta6 = 0

    kappa = 40

    l = 1

    print(ComputeRuntimeAndAlpha(n, t, circuit_size, 100, 20, d,
                           beta1, beta2, beta3, beta4, beta5, beta6, kappa))

    for l in (1,5,10,20):
        ComputeRuntimesAndAlphasForK(n, t, circuit_size, l, d,
                                 beta1, beta2, beta3, beta4, beta5, beta6, kappa)

