#!/usr/bin/env python3
"""
    This code generates hyper-surface for simple geometry to test particle
    samplers
"""

import numpy as np

def Bjorken_hydro(tau_fo, T_fo):
    """This function outputs a hyper-surface from a 0+1D Bjorken hydro"""
    R_fo = 10.      # fm
    dx   = 0.1      # fm
    dtau = 0.1      # fm
    deta = 1.0
    volume = tau_fo*dtau*dx*dx*deta

    nx  = 2.*R_fo/dx + 1
    x   = np.linspace(-R_fo, R_fo, nx)
    y   = np.linspace(-R_fo, R_fo, nx)
    R   = np.sqrt(x**2. + y**2.)
    idx = R <= R_fo

    xfo       = x[idx]
    yfo       = y[idx]
    fo_length = xfo.shape
    taufo     = tau_fo*np.ones(fo_length)
    eta       = np.zeros(fo_length)

    dsigma_0 = volume*np.ones(fo_length)
    dsigma_1 = np.zeros(fo_length)
    dsigma_2 = np.zeros(fo_length)
    dsigma_3 = np.zeros(fo_length)

    utau_fo = np.ones(fo_length)
    ux_fo   = np.zeros(fo_length)
    uy_fo   = np.zeros(fo_length)
    ueta_fo = np.zeros(fo_length)

    e_fo = T_fo**4.*np.ones(fo_length)
    Tfo  = T_fo*np.ones(fo_length)
    muBfo = np.zeros(fo_length)
    muSfo = np.zeros(fo_length)
    muQfo = np.zeros(fo_length)
    eplusp_over_T_FO = T_fo**3.*np.ones(fo_length)

    Wtautau = np.zeros(fo_length)
    Wtaux   = np.zeros(fo_length)
    Wtauy   = np.zeros(fo_length)
    Wtaueta = np.zeros(fo_length)
    Wxx     = np.zeros(fo_length)
    Wxy     = np.zeros(fo_length)
    Wxeta   = np.zeros(fo_length)
    Wyy     = np.zeros(fo_length)
    Wyeta   = np.zeros(fo_length)
    Wyeta   = np.zeros(fo_length)
    Wetaeta = np.zeros(fo_length)

    pi_b  = np.zeros(fo_length)
    rho_b = np.zeros(fo_length)
    qtau  = np.zeros(fo_length)
    qx    = np.zeros(fo_length)
    qy    = np.zeros(fo_length)
    qeta  = np.zeros(fo_length)

    output = np.array([taufo, xfo, yfo, eta,
                       dsigma_0, dsigma_1, dsigma_2, dsigma_3,
                       utau_fo, ux_fo, uy_fo, ueta_fo, e_fo, Tfo,
                       muBfo, muSfo, muQfo, eplusp_over_T_FO,
                       Wtautau, Wtaux, Wtauy, Wtaueta, Wxx, Wxy, Wxeta,
                       Wyy, Wyeta, Wetaeta, pi_b, rho_b, qtau, qx, qy, qeta])
    print(output)
    output.transpose().astype('float').tofile("frzout.dat")

if __name__ == "__main__":
    Bjorken_hydro(10., 0.15)
