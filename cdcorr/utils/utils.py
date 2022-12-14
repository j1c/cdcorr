import numpy as np

def example1(n):
    mean = np.array([0, 0, 0])
    cov = np.array([
        [1, 0.36, 0.6],
        [0.36, 1, 0.6],
        [0.6, 0.6, 1]]
    )
    
    X = np.random.multivariate_normal(mean, cov, size=n)
    Y = np.random.multivariate_normal(mean, cov, size=n)
    Z = np.random.multivariate_normal(mean, cov, size=n)
    
    return X, Y, Z


def compute_kernel(X, sigma):
    """
    X : (n_samples, n_dim)
    sigma : float
    """

    d = X.shape[1]
    denom = np.power(2 * np.pi, d / 2.0) * np.power(sigma, d / 2)
    constant = 1 / denom

    kern = pairwise_kernels(X, metric="rbf", gamma=1 / (sigma * 2)) * constant

    return kern