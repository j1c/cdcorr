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