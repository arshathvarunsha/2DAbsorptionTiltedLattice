# jit
import time
import matplotlib.pylab as plt
import numpy as np
from scipy.fft import fft, ifft
# print("Using splitop_tilted.py")
# input()
# from numba import jit, objmode
class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

# Generating the initial state of the electronic part:
def initElectronic(parameters, dat, initState = 0):
    initype = parameters.initype
    Nsite = parameters.Nsite
    Nlayer = parameters.Nlayer
    Npht = parameters.N_Pht
    if initype == 4:
        c = np.zeros(Nsite * Nlayer + Npht) + 0.j
        c[initState] = 1/np.sqrt(2)      
        c[initState + 1] = 1/np.sqrt(2) *np.exp(1j*parameters.θr)
        fact = parameters.Nsite*parameters.Nlayer
        print(parameters.ωk[initState-fact]*27.2114, "ωk")
        print(parameters.ωk[initState - fact + 1]*27.2114, "ωk")
        print(parameters.ωk[initState-fact + 2*parameters.δk]*27.2114, "ωk")
        print(parameters.ωk[initState - fact + 1+ 2*parameters.δk]*27.2114, "ωk")
        print(parameters.ωk[0]*27.2114, "ωk0")
    if max(abs(c))==0 or max(abs(c))==float('nan'):
        c[initState] = 1
        exit()
    return c


def Uγ(Ψ, R, param, dt):
    ϵ0 = param.ϵ0  - 1j * param.Γ/2
    Ψ[:param.Nsite * param.Nlayer] = Ψ[:param.Nsite * param.Nlayer] * np.exp(-1j * (ϵ0 + param.γ * R ) * dt)
    return Ψ

def Uk(Ψ, param, dt):
    # Matter
    matter = param.Nsite * param.Nlayer
    Ψ_layers = Ψ[:matter] * 1.0
    Ψ_layers = Ψ_layers.reshape(param.Nlayer, param.Nsite)
    Ψfft = np.fft.fft(Ψ_layers, axis=1, norm='ortho') 
    Ψm = Ψfft.flatten()
    k = np.fft.fftfreq(param.Nsite, d= 1/param.Nsite)
    Ek = (-2) * param.τ * np.cos((2 * np.pi * k)/param.Nsite)
    Ek = np.tile(Ek, param.Nlayer)
    ψ_m = Ψm * np.exp(-1j * Ek * dt)
    Ψ[matter:] *= np.exp(-1j * param.ωk * dt)
    #for i in range(len(Ψ[matter:])):   
    ψ_pi = Ψ[matter:] *1.0
    Ψnl = np.reshape(ψ_m, (param.Nlayer, param.Nsite))
    Ψr  = np.fft.ifft(Ψnl, norm='ortho', axis=1)
    Ψr = Ψr.flatten()
    ψ_pi = np.fft.ifft(ψ_pi * 1.0, norm='ortho')
    Ψr = np.concatenate((Ψr * 1.0, ψ_pi * 1.0))
    Ψr = Ψr.flatten()
    return Ψr

def Uc(Ψ, param, dt):
    Nx = param.Nsite
    Np = param.N_Pht
    L = param.Nlayer
    Srs =  np.sin(param.kz *  param.Z_nm)
    Nn = (np.sum((Srs**2), axis=1))
    eta = param.g * (np.sqrt(Nn))
    Srs = np.einsum("nl, n -> nl", Srs, 1/np.sqrt(Nn))
    Ψexi = np.reshape(Ψ[:param.Nsite*param.Nlayer] *1.0, (param.Nsite, param.Nlayer))
    Ψex =  np.sum((Srs*Ψexi), axis=1)
    Ψex =  Ψex.flatten()
    Ψph =  Ψ[param.Nsite*param.Nlayer:]  * 1.0
    Ψph =  Ψph.flatten()

    Ep  =     eta
    Em  = -   eta
    Ψp =  (np.exp(-1j * Ep * dt) * (1/np.sqrt(2)*Ψex + 1/np.sqrt(2)*Ψph))
    Ψm =  (np.exp(-1j * Em * dt) * (1/np.sqrt(2)*Ψex - 1/np.sqrt(2)*Ψph))
    Ψex_ = 1/np.sqrt(2) * (Ψp + Ψm)
    Ψph_ = 1/np.sqrt(2) * (Ψp - Ψm)


    δΨex = Ψex_ - Ψex
    Ψex_nl = np.einsum("nl,n->nl", Srs, δΨex)
    δΨc = Ψex_nl.flatten()
    Ψf = np.zeros(Nx*L+Np)  + 0j
    Ψf[0:Nx*L] = Ψ[0:Nx*L] * 1.0 + δΨc
    Ψph_ = np.fft.fft(Ψph_ * 1.0, norm='ortho')   
    Ψf[Nx*L:] = Ψph_
    Ψf = Ψf.flatten()
    return Ψf

def propagateCi(parameters, dat, ci, dt):
    ci = ci * 1.0 + 0j  
    ct_trot = Uγ(ci, dat.R, parameters, dt/2) * 1.0    
    ct_trot = Uk(ct_trot, parameters, dt) * 1.0
    ct_trot = Uc(ct_trot, parameters, dt) * 1.0
    ct_trot = Uγ(ct_trot, dat.R, parameters, dt/2) * 1.0
    return ct_trot




#@jit(nopython=False)
def Force(γ, dH0, ci):
    F = -dH0  
    ndof = len(dH0)
    F -= γ * (ci.conj() * ci).real[:ndof] 
    return F

def VelVer(dat) : 
    par =  dat.param
    v = dat.P/par.M
    F1 = dat.F1 
    # electronic wavefunction
    ci = dat.ci * 1.0
    dat.ck_1 =  dat.ck_1 * 0.0
    dat.ck_2 =  dat.ck_2 * 0.0
    dat.ck_δ =  dat.ck_δ * 0.0
    dat.ck_δ2 =  dat.ck_δ2 * 0.0
    # ======= Nuclear Block ==================================
    EStep = int(par.dtN/par.dtE)
    dtE = par.dtN/EStep
    # half electronic evolution
    for t in range(int(np.floor(EStep/2))):
        dat.ck_1[t] = ci[par.initState]
        dat.ck_2[t] = ci[par.initState + 1]
        dat.ck_δ[t] = ci[par.initState + par.δk * 2]
        dat.ck_δ2[t] = ci[par.initState + 1 + par.δk * 2]
        ci = propagateCi(par, dat, ci, dtE)         
    #ci /= np.sum(ci.conjugate()*ci) 
    dat.ci = ci * 1.0 
    # ======= Nuclear Block ==================================
    dat.R += v * par.dtN + 0.5 * F1 * par.dtN ** 2 / par.M  
    #------ Do QM ----------------
    #dat.Hij  = par.Hel() + 0j  
    dat.Hij  = par.HelR(dat.R)  
    dat.dH0  = par.dHel0(dat.R)
    #-----------------------------
    F2 = Force(par.γ, dat.dH0, dat.ci) # force at t2
    v += 0.5 * (F1 + F2) * par.dtN / par.M
    dat.F1 = F2
    dat.P = v * par.M
    # ======================================================
    # half electronic evolution
    for t in range(int(np.ceil(EStep/2))):
        dat.ck_1[t + int(np.floor(EStep/2))] = ci[par.initState]
        dat.ck_2[t + int(np.floor(EStep/2))] = ci[par.initState + 1]
        dat.ck_δ[t + int(np.floor(EStep/2))] = ci[par.initState + par.δk * 2]
        dat.ck_δ2[t + int(np.floor(EStep/2))] = ci[par.initState + 1 + par.δk * 2]
        ci = propagateCi(par, dat, ci, dtE) 
    dat.ci = ci * 1.0 

    return dat


def pop(dat):
    ci =  dat.ci
    return np.outer(ci.conjugate(),ci)

def runTraj(parameters):
    #------- Seed --------------------
    try:
        np.random.seed(parameters.SEED)
    except:
        pass
    #------------------------------------
    ## Parameters -------------
    dat = Bunch(param =  parameters )
    NSteps = parameters.NSteps
    NTraj = parameters.NTraj
    Nsite = parameters.Nsite
    Nlayer = parameters.Nlayer
    initState = parameters.initState # intial state
    nskip = parameters.nskip
    mask = parameters.mask
    #---------------------------
    if NSteps%nskip == 0:
        pl = 0
    else :
        pl = 1    
    Esteps = int(parameters.dtN/parameters.dtE)
    ## Initializing
    rho_ensemble = np.zeros(((NSteps//nskip + pl)*Esteps), dtype=complex)
    rho_ensemble1 = np.zeros(((NSteps//nskip + pl)*Esteps), dtype=complex)
    rho_ensemble2 = np.zeros(((NSteps//nskip + pl)*Esteps), dtype=complex)
    rho_ensembleδ = np.zeros(((NSteps//nskip + pl)*Esteps), dtype=complex)
    rho_ensembleδ2 = np.zeros(((NSteps//nskip + pl)*Esteps), dtype=complex)
    #---------------------------
    # Parameters
    dat.Pht_ks = np.arange(Nsite)[mask]
    tinit0 = time.time()
    # Call function to initialize mapping variables
    print("Initial state:", time.time() - tinit0)
    for itraj in range(NTraj):
        seednum = parameters.SEED + itraj + int(time.time() * 0.00001)
        print("Seednum: ", seednum)
        # Trajectory data
        dat.R, dat.P = parameters.initR(seednum)
        dat.ci = initElectronic(parameters, dat, initState) # np.array([0,1])
        # set propagator
        vv  = VelVer
        #----- Initial QM --------
        dat.Hij  = parameters.HelR(dat.R)  
        dat.dH0  = parameters.dHel0(dat.R)
        dat.F1 = Force(parameters.γ, dat.dH0, dat.ci) # Initial Force
        #----------------------------
        iskip = 0 # please modify
        t0 = time.time()
        dat.ck_1 = np.zeros((Esteps)) + 0j
        dat.ck_2 = np.zeros((Esteps)) + 0j
        dat.ck_δ = np.zeros((Esteps)) + 0j
        dat.ck_δ2 = np.zeros((Esteps)) + 0j
        for i in range(NSteps): # One trajectory
            dat = vv(dat)
            #------- ESTIMATORS-------------------------------------
            if (i % nskip == 0):
                rho_ensemble1[iskip * Esteps : (iskip + 1) * Esteps] += dat.ck_1[:] * 1.0
                rho_ensemble2[iskip * Esteps : (iskip + 1) * Esteps] += dat.ck_2[:] * 1.0
                rho_ensembleδ[iskip * Esteps : (iskip + 1) * Esteps] += dat.ck_δ[:] * 1.0
                rho_ensembleδ2[iskip * Esteps : (iskip + 1) * Esteps] += dat.ck_δ2[:] * 1.0
                iskip += 1
            #-------------------------------------------------------        
        time_taken = time.time()-t0
        print(f"Time taken: {time_taken} seconds")


    return rho_ensemble1, rho_ensemble2, rho_ensembleδ, rho_ensembleδ2