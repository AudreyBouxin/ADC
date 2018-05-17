# -*- coding: utf-8 -*-
"""
Created on Thu May 17 13:37:34 2018

@author: audrey.bouxin
"""

import numpy as np
import matplotlib.pyplot as plt


zenith_angle = 45   #[°]
lambda_wave = 1e-9*np.linspace(300,900,100) #np.array([300e-9,900e-9])  #[m]

#constantes
DEG2RAD = np.pi/180

def fRefractionAtm():
    """
    ## Fonction qui calcule la réfraction de l'atmosphère pour des longueurs d'onde données
    #
    # INPUTS :
    #   - lambda       [m] vecteur de longueurs d'onde étudiées
    #   - zenith_angle [°] l'angle au zénith sous lequel on observe
    # OUTPUT :
    #   - Ratm       [rad] la dispersion angulaire de l'atmosphère
    """
    #constantes:
    mmHg2Pa = 133.322365;
    zenith_angle_rad = DEG2RAD*zenith_angle;
    lambda_wave_um = lambda_wave*1e6;

    # Initialisation du problème
    #Constants involved in the standard phase and group refractivities of dry air (Ciddar 1996)
    k0 = 238.0182;	#[um^-2]
    k1 = 5792105;   #[um^-2]
    k2 = 57.362;    #[um^-2]
    k3 = 167917;    #[um^-2]
    
    #Constants involved in the standard phase and group refractivities of water vapor (Ciddar 1996)
    w0 = 295.235;   #[um^-2]
    w1 = 2.6422;    #[um^-2]
    w2 = -0.032380;	#[um^-4]
    w3 = 0.004028;  #[um^-6]
    
    #les équations d calcul de l'indice de réfraction de l'atmosphère pour le
    #visible et l'infra rouge proche sont alors
    sigma = 2*np.pi/lambda_wave_um;  #[m^-1] le nombre d'onde
    
    nas = 1+1e-8*(k1/(k0-sigma**2)+k3/(k2-sigma**2));  #refractive index of standard air at 15°C, 101325Pa,0#humidity, 450ppm of CO2
    
    #For us it is not usefull to have the following equation because we are in
    #a standard case so naxs=nas, but I let it in order to remember the sequence
    #of equations
    xc = 450;                                   #[ppm] number of ppm of CO2 (standard) (Ciddar 1996)
    naxs = 1+(nas-1)*(1+0.534e-6*(xc-450));    #refractive index of standard air at 15°C, 101325Pa,0#humidity, xc ppm of CO2
    
    nws = 1+1.022e-8*(w0+w1*sigma**2 +w2*sigma**4+w3*sigma**6); #refractive index of water vapor at 20°C, 1335Pa
    
    #Paramètre généraux
    Height = 3170;                  #[m] altitude Karakaya
    T  = -15;                       #[°C] température
    T_K = T + 273.15;               #[K] température
    T0  = 15+273.15;                #[K] température standard (Robo-AO)
    DT  = 0.0065;                   #[K] gradient vertical de température 0.65K pour 100m (Wikipedia)
    P0 = 1.01325e5;                 #[Pa] pression normal à l'altitude 0
    P  = P0*(1-DT*Height/T_K);      #[Pa] pression (Wikipedia)
    RH = 0.8;                       #[-] humidité relative
    
    xCO2 = 0.0004;                              #(Davis 1992)
    Ma = (28.9635+12.011*(xCO2-0.0004))*1e-3;   #[kg/mol] (Davis 1992)
    Mv = 18.01528e-3;                           #[kg/mol] (Google)
    R = 8.314510;                               #[J/mol/K] (Davis 1992)
    
    #Coefficients pour pression de vapeur saturante Psv (Davis 1992)
    A = 1.2378847e-5;   #[K^-2]
    B = -1.9121316e-2;  #[K^-1]
    C = 33.93711047;    #[-]
    D = -6.4341645e3;   #[K]
    
    #Coefficients pour facteur d'augmentation f (Davis 1992)
    alpha = 1.00062;    #[-]
    beta  = 3.14e-8;    #[Pa^-1]
    gamma = 5.6e-7;     #[K^-2]
    
    #Coefficients pour facteur de compressibilité Z (Davis 1992)
    a0 = 1.58123e-6; #[K*Pa^-1]
    a1 = -2.9331e-8; #[Pa^-1]
    a2 = 1.1043e-10; #[(K*Pa)^-1]
    b0 = 5.707e-6;   #[K*Pa^-1]
    b1 = -2.051e-8;  #[Pa^-1]
    c0 = 1.9898e-4;  #[K*Pa^-1]
    c1 = -2.376e-6;  #[Pa^-1]
    d  = 1.83e-11;   #[K^2Pa^-2]
    e  = -0.765e-8;  #[K^2Pa^-2]
    
    Psv = np.exp(A*T_K**2+B*T_K+C+D/T_K);   #la pression de vapeur d'eau saturante dans l'air humide (Davis 1992)
    f = alpha+beta*P+gamma*T_K**2;       #facteur d'augmentation (Davis 1992)
    xv = RH*f*Psv/P;                    #fraction molaire de la vapeur d'eau
    #Z = 1-P/T_K*(a0+a1*T_K+a2*T_K**2+(b0+b1*T_K)*xv +(c0+c1*T_K)*xv**2)
    #+(d+e*xv**2)*(P/T_K)**2;   #compressibilité erreur sur les T qui doivent
    #être en °C pour certains
    Z = 1-P/T_K*(a0+a1*T+a2*T**2+(b0+b1*T)*xv +(c0+c1*T)*xv**2) +(d+e*xv**2)*(P/T_K)**2;   #compressibilité
    ro_a = P*Ma/(Z*R*T_K)*(1-xv*(1-Mv/Ma)); #[kg/m^3] masse volumique de l'air humide
    
    Rgas = 287.05;              #[J/kg/K] https://www.brisbanehotairballooning.com.au/calculate-air-density/
    ro_axs = P0/(Rgas*T0);      #[kg/m^3] masse volumique de l'air sec dans les conditions standards
    
    # Définition des variables de pression
    Psat = np.exp(46.784-6435/(T_K)-3.868*np.log(T_K));	#[mmHg] valide entre -50°C et 200°C
    Pw = Psat*RH*mmHg2Pa;                           #[mb] pression partielle de vapeur d'eau
    Psat0 = np.exp(46.784-6435/(T0)-3.868*np.log(T0));	#[mmHg] valide entre -50°C et 200°C
    Pw0 = Psat0*RH*mmHg2Pa;                         #[Pa] pression partielle de vapeur d'eau
    #Ps = P-Pw;                                      #[Pa] pression partielle d'air sec
    ro_w = Pw*Mv/(R*T_K);                           #[kg/m^3] masse volumique de la vapeur d'eau
    ro_ws = Pw0*Mv/(R*T0);                          #[kg/m^3] masse volumique de la vapeur d'eau dans les conditions standards
    
    # Ciddor equations for the refractive index of the atmosphere
    n = 1+(ro_a/ro_axs)*(naxs-1)+(ro_w/ro_ws)*(nws-1);
    
    
    # La réfraction de l'atmosphère (Robo-AO)
    kappa = 1;  #constante cf(Robo-AO)
    beta = 0.001254*(T_K/273.15);   #cf(Robo-AO)
    #on peut négliger le term en ^3 pour des angles au zénith < 65°
    Ratm = kappa*(n-1)*(1-beta)*np.tan(zenith_angle_rad)-kappa*(n-1)*(beta-(n-1)/2)*(np.tan(zenith_angle_rad))**3
    return Ratm



def plotRatm():    
#    from fRefractionAtmosphere import fRefractionAtm
    Ratm= fRefractionAtm()
    plt.figure()
    plt.plot(lambda_wave*1e9,Ratm)  