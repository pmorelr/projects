import numpy as np

''''Comentarios Gerais'''

# os outputs de todas as funcoes retornam uma lista 'return'
# essa lista contem 4 tipos de resultados diferentes:
# [erosion rate , volume loss , relative thickness loss, annual thickness loss]
# em geral, as funcoes referentes aos modelos de erosao  classicos' possuirao os dois primeiros resultados
# os dois ultimos sao mais comuns aos modelos de erosao por geometria, da DNV
# a fim de representar a ausencia de um tipo de resultado, usa-se o 0

''''Variaveis'''
# vp    -> particle velocity        (m/s)
# alfa  -> impact angle             (degrees)
# dp    -> particle diameter        (m)
# rot   -> target wall density      (kg/cm3)
# mp    -> particle mass            (kg)
# Hv    -> vickers hardness         ()
# rop   -> particle density         (kg/cm3)
# D     -> pipe diameter            (m)
# D1/D2 -> pipe diameter (reducer)  (m)
# mum   -> mixture fluid viscosity  (kg/ms)
# rom   -> mixture fluid density    (kg/m3)
# R     -> pipe bend radius         (pipe diameters)
# BH    -> brinell hardness         ()
# P     -> plastic flow stress      (Pa)
# G     -> geometry factor          (1, 2, 3. ver documento da dnv)
# Fs    -> sharpness factor         (0.2 -> fully rounded / 0.5 -> semi-rounded / 1.0 -> sharp)



''' Impinging Jet Models '''

# Ahlert
def ahlert(vp, alfa_degree, rot, mp, H, Fs):
    alfa = np.radians(alfa_degree)

    g_alfa = 5.40*alfa + -10.11*alfa**2 + 10.93*alfa**3 + -6.33*alfa**4 + 1.42*alfa**5
    er = 2.17e-7*H**(-0.59)*Fs*vp**(2.41)*g_alfa

    volume_loss = er*mp/rot

    return [er,volume_loss, 0, 0]

# Arabnejad
def arabnejad(vp, alfa_degree, rot, mp, Fs):
    alfa = np.radians(alfa_degree)

    # df -> deformation wear factor
    df = 3.5e11
    # C -> cutting constant
    C = 0.015

    # Values for stainless steel 2205
    # vpt -> threshold velocity (velocidade minima para que a deformacao seja considerada)
    vpt = 2.3
    # k -> particle shape/material deformation constant
    k = 0.4
    c1 = 3.92e-8
    c2 = 2.30e-8

    listlow = c1*Fs*vp**2.41*np.sin(alfa)*(2*k*np.cos(alfa)-np.sin(alfa))/(2*k**2)
    listhigh = c1*Fs*vp**2.41*np.cos(alfa)**2/2

    erc = np.where(alfa<np.arctan(k),listlow,listhigh)

    erd = c2*Fs*(max(vp*np.sin(alfa) - vpt, 0))**2

    er = erc + erd

    volume_loss = er*mp/rot

    return [er,volume_loss, 0, 0]

# DNV
def dnv(vp, alfa_degree, rot, mp):
    alfa = np.radians(alfa_degree)

    # empirical constants
    n=2.6
    K=2.0e-9

    H_alfa = 0.6*(np.sin(alfa)+7.2*(np.sin(alfa)- np.sin(alfa)**2))**0.6*(1-np.exp(-20*alfa))
    er = K*vp**n*H_alfa
    volume_loss = er*mp/rot
    return [er,volume_loss, 0, 0]

# Finnie
def finnie(vp, alfa_degree,rot, mp, P):
    alfa = np.radians(alfa_degree)

    # k -> ratio of vertical to horizontal force component on particle
    k = 2.0
    # psi -> ratio of depth of contact to depth of cut
    psi = 1

    if np.tan(alfa)<=k/6.0:
        er = 0.5*((mp*vp**2)/(P*psi) * (np.sin(2*alfa) - 6*np.sin(alfa)**2/k))
    else:
        er = 0.5*((mp*vp**2)/(P*psi) * ((k*np.cos(alfa)**2)/6))

    result = er/mp*rot
    volume_loss = er*mp/rot

    return [result,volume_loss,0,0]

# Haugen
def haugen(vp, alfa_degree, rot, mp):
    alfa = np.radians(alfa_degree)

    # empirical constants
    k = 2e-9
    n = 2.6

    f_alfa = 9.370*alfa - 42.295*alfa**2 + 110.804*alfa**3 - 175.804*alfa**4 + 170.137*alfa**5 - 98.298*alfa**6 + 31.211*alfa**7 - 4.170*alfa**8

    er = mp*k*f_alfa*vp**n
    volume_loss = er*mp/rot

    return [er,volume_loss, 0, 0]

# Huang
def huang(vp, alfa_degree, dp, rot, mp, rop):
    alfa = np.radians(alfa_degree)

    # C -> cutting constant
    C= 7.5e-3
    # D -> deformation constant
    D= 0.082

    er1 =  C*mp*rop**0.15*(vp*np.sin(alfa))**2.3

    er2 = D*mp*rop**0.1875*dp**(0.5)*vp**2.375*np.cos(alfa)**2*np.sin(alfa)**0.375
    total = er1 + er2
    er = total*rot/mp *1e-9
    volume_loss = er*mp/rot
    return [er,volume_loss, 0, 0]

# Oka
def oka(vp, alfa_degree, dp, rot, mp, Hv):
    alfa = np.radians(alfa_degree)
    Hv_GPa = 0.009807 * Hv

    # empirical constants
    k3=0.19
    k1= -0.12
    Kp=65
    # reference particle diameter
    dl= 326e-6
    # reference particle speed
    vl=104

    n1 = 0.71*(Hv_GPa)**0.14
    n2 = 2.4*(Hv_GPa)**(-0.94)

    k2 = 2.3*(Hv_GPa)**0.038

    e90 = Kp*(Hv_GPa)**(k1)*(vp/vl)**k2*(dp/dl)**k3
    falfa = (np.sin(alfa))**n1*(1+Hv*(1-np.sin(alfa)))**n2
    ev_alfa = falfa*e90
    er = 1e-9*rot*ev_alfa

    volume_loss = er*mp/rot
    return [er,volume_loss, 0, 0]

''' DNV Models '''

# Smooth Pipe
def smooth(vp, mp, D):
    rel_thick_loss = 8e-10*vp**2.6*D**(-2)
    thick_loss = rel_thick_loss*mp*1e-6
    return [0, 0,rel_thick_loss,thick_loss]

# Welded Joint
def welded(vp, alfa_degree, dp, rot, mp, D, rom):
    alfa = np.radians(alfa_degree)

    # empirical constants
    n=2.6
    K=2.0e-9

    c_factor = 1e6*dp/(30*rom**0.5)
    c_factor = min(c_factor, 1.0)
    H_alfa = 0.6*(np.sin(alfa)+7.2*(np.sin(alfa)- np.sin(alfa)**2))**0.6*(1-np.exp(-20*alfa))

    A_pipe = np.pi*D**2/4
    rel_thick_loss = K*vp**n*H_alfa*np.sin(alfa)/(rot*A_pipe)*c_factor*1e6
    thick_loss = rel_thick_loss*mp*1e-6

    return [0, 0,rel_thick_loss,thick_loss]

# Pipe Bend
def bend(vp,dp, rot, mp, rop, D, mum, rom, R, GF):
    alfa = np.arctan(1/(np.sqrt(2*R)))
    H_alfa = 0.6*(np.sin(alfa)+7.2*(np.sin(alfa)- np.sin(alfa)**2))**0.6*(1-np.exp(-20*alfa))
    A = rom**2*np.tan(alfa)*vp*D/(rop*mum)
    B = rop/rom
    gamac = 1/(B*(1.88*np.log(A)-6.04))
    gamac = min(gamac, 0.1)

    gama = dp/D

    if gama < gamac:
        G = gama/gamac
    else:
        G = 1

    area = np.pi*D**2/(4*np.sin(alfa))

    C1 = 2.5

    K=2.0e-9
    n=2.6
    rel_thick_loss = K*vp**n*H_alfa*G*C1*GF*1e6/(rot*area)
    thick_loss = rel_thick_loss*mp*1e-6

    return [0, 0,rel_thick_loss, thick_loss]

# Blinded Tee
def blinded(vp, dp, rot, mp, rop, D, mum, rom, GF):

    gama = dp/D
    beta = rop/rom
    Re_d = vp*D*rom/mum
    if beta>40:
        gamac = 0.14/beta
        if gama < gamac:
            c = 19/(np.log(Re_d))
        else:
            c = 0
        C1 = 3/(beta**0.3)
    else:
        b = (np.log(Re_d/10000 + 1) + 1)**(-0.6) -1.2
        gamac = 0.0035*(beta/40)**b
        if gama < gamac:
            c = 19/(np.log(Re_d))
        else:
            c = -0.3*(1 - 1.01**(40 - beta))
        C1 = 1.0

    G = (gama/gamac)**c

    At = np.pi*D**2/4

    K=2.0e-9
    n=2.6
    rel_thick_loss = K*vp**n/(rot*At)*G*C1*GF*1e6
    thick_loss = rel_thick_loss*mp*1e-6
    return [0, 0,rel_thick_loss, thick_loss]


# Reducer
# obs: vp1 se refere a velocidade da particula no tubo com diametro D1.
def reducer(vp1, alfa_degree, dp, rot, mp, rom, D1, D2, GF):
    alfa = np.radians(alfa_degree)

    H_alfa = 0.6*(np.sin(alfa)+7.2*(np.sin(alfa)- np.sin(alfa)**2))**0.6*(1-np.exp(-20*alfa))
    At = np.pi*(D1**2 - D2**2)/(4*np.sin(alfa))
    Aratio = 1 - (D2/D1)**2
    vp = vp1*(D1/D2)**2
    fator = 1e6*dp/(30*rom**0.5)
    C2 = min(fator, 1.0)

    K=2.0e-9
    n=2.6
    rel_thick_loss =  H_alfa*K*vp**n/(rot*At)*Aratio*C2*GF*1e6
    thick_loss = rel_thick_loss*mp*1e-6

    return [0, 0,rel_thick_loss, thick_loss]
