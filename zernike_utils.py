"""
Utility functions for Zernike polynomial generation and analysis.

This module provides functions to compute Zernike polynomials using
Noll indices, convert mode amplitudes to normalized units, and evaluate
modes over polar coordinates. It also includes a helper function to
compute the root-mean-square (RMS) of data arrays.

Intended Use
------------
- Generating target Zernike surfaces for deformable mirror FEM simulations.
- Evaluating actuator influence function projections against Zernike modes.
- Quantifying residual errors in optical surface fitting.

Functions
---------
- zernike_amplitude(weight) : Scale a Zernike mode amplitude (nm) for FEM.
- zernike_radial(n, m, r) : Compute the radial polynomial R_n^m(r).
- noll_to_nm(j) : Convert a Noll index to radial (n) and azimuthal (m) orders.
- zernike_noll(R, Phi, j, weight) : Evaluate a Zernike mode over (R, Phi).
- rms(data) : Compute root-mean-square of an array.


"""



import numpy as np
import math





def zernike_amplitude(weight):             

    norm=(1/4000)*weight        # weighted Zernike amplitud (250 nm) 
     
    return norm
    
    
    
    
    
    

def zernike_radial(n, m, r):


    """
    Compute Zernike radial polynomial R_n^m(r).

    Parameters
    ----------
    n : int
        Radial order
    m : int
        Azimuthal order (>= 0)
    r : ndarray
        Normalized radial coordinate (0 <= r <= 1)
    """

    R = np.zeros(r.shape)

    for i in range(0, int((n - m) / 2) + 1):

        R += np.array(r**(n - 2 * i) * (((-1)**(i)) *
                         np.math.factorial(n - i)) /
                         (np.math.factorial(i) *
                          np.math.factorial(int(0.5 * (n + m) - i)) *
                          np.math.factorial(int(0.5 * (n - m) - i))),
                         dtype='float')
    return R




def noll_to_nm(j):
    """
    Convert Noll index to (n, m).

    Parameters
    ----------
    j : int
        Noll index (j >= 1)

    Returns
    -------
    n : int
        Radial order
    m : int
        Azimuthal order (can be negative)
    """

    if j < 1:
        raise ValueError("j must be >= 1")

    # 1. Radial order
    n = int(math.floor((math.sqrt(8*j - 7) - 1) / 2))

    # 2. Index within order
    k = j - n*(n+1)//2 - 1

    # 3. Allowed |m|
    m_abs = list(range(n % 2, n + 1, 2))

    # 4. Sign ordering (Noll)
    if n % 4 in (0, 1):
        signs = [1, -1]
    else:
        signs = [-1, 1]

    m_list = []
    for ma in m_abs:
        if ma == 0:
            m_list.append(0)
        else:
            for s in signs:
                m_list.append(s * ma)

    return n, m_list[k]





def zernike_noll(R, Phi, j, weight):
    """
    Compute Zernike polynomial using Noll index.

    Parameters
    ----------
    R : ndarray
        Normalized radial coordinate
    Phi : ndarray
        Azimuthal coordinate (radians)
    j : int
        Noll index
    weight : float
        Mode amplitude scaling

    Returns
    -------
    Z : ndarray
        Zernike mode evaluated at (R, Phi)
    """
    n, m = noll_to_nm(j)
    norm = zernike_amplitude(weight)

    m_abs = abs(m)
    radial = zernike_radial(n, m_abs, R)

    if m == 0:
        Z = radial
    elif m > 0:
        Z = radial * np.cos(m_abs * Phi)
    else:
        Z = radial * np.sin(m_abs * Phi)

    return norm * Z




      
  
def rms(data):
	return np.sqrt(np.mean(np.square(data)))


