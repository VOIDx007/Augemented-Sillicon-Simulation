import numpy as np
clear
## Constant definition cell
hbar = 1.06e-34
q = 1.6e-19
epsil = 10 * 8.85e-12
kTH = 8.314 / 6.023e+23 * 330 / q

kTC = 8.314 / 6.023e+23 * 300 / q

kT = 0.5 * (kTH + kTC)

m = 0.07 * 9.1e-31
n0 = 2 * m * kT * q / (2 * np.pi * (hbar ** 2))
IE = (q * q) / (2 * np.pi * hbar)
n0H = 2 * m * kTH * q / (2 * np.pi * (hbar ** 2))
n0C = 2 * m * kTC * q / (2 * np.pi * (hbar ** 2))
tic
a = 2.5e-10

t = (hbar ** 2) / (2 * m * (a ** 2) * q)
beta = q * a * a / epsil
## Device definition
#Device is divided into 3 parts: Heavily doped contacts (source and drain),
# lightly doped active region and lightly doped spacer (between contacts
# and active region). The active region consists of the two barriers and
# well. The active region together with the spacers is referred to as the
# channel.

#The "quantum" simulation is only performed in the channel region, since
#confinement only exists here. In the contact regions a semiclassical
#equation is used for electron population, as described later. Quantum
#simulation for the whole device could be used, but it would lead to a much
#larger computational time.

#In principle the effective mass and permittivity should be different for
#barrier (AlAs) and well (GaAs). however since the electron concentration
#in the barrier is really small, and the permittivities of GaAs and AlAs
#are similar, both effective mass and permittivity are taken to be constant
#throughout the device. This should, however, be corrected in the long
#run, and would require manual summing over all lateral modes (with
#constant effective mass throughout, this summing results in the F2D
#functions being used for electron population instead of the fermi
#function). I have done this for a simpler problem (self-consistent
#Schrodinger-Poission solver); please contact me in case manual sum over
#modes proves problematic.
#Only a small effect is expected on calculations.
Ls = 3e-08

Lbl = 4e-09

Lbr = 4e-09

Lg = 7e-09

Lbuff = 1e-09

Ns = np.round(Ls / a)
Nbl = np.round(Lbl / a)
Nbr = np.round(Lbr / a)
Ng = np.round(Lg / a)
Nl = np.round(2 * Lbuff / a)

Nc = Nl + Nbl + Nbr + Ng

Np = Ns + Nc + Ns
XX = a * 1000000000.0 * np.array([np.arange(1,Np+1,1)])
mu = 0
N_cH = 2 * ((n0H / 2) ** 1.5)

N_cC = 2 * ((n0C / 2) ** 1.5)

N_c = 2 * ((n0 / 2) ** 1.5)

Nd = 1e+24

Nd = Nd * np.array([[np.ones((Ns,1))],[0.001 * np.ones((Nc,1))],[np.ones((Ns,1))]])

D2 = - (2 * diag(np.ones((1,Np)))) + (diag(np.ones((1,Np - 1)),1)) + (diag(np.ones((1,Np - 1)),- 1))

D2[1,1] = - 1
D2[Np,Np] = - 1

T = (2 * t * diag(np.ones((1,Np)))) - (t * diag(np.ones((1,Np - 1)),1)) - (t * diag(np.ones((1,Np - 1)),- 1))
Ec = np.array([[np.zeros((Ns,1))],[np.zeros((Nl / 2,1))],[0.1 * np.ones((Nbl,1))],[np.zeros((Ng,1))],[0.1 * np.ones((Nbr,1))],[np.zeros((Nl / 2,1))],[np.zeros((Ns,1))]])
T = T + diag(Ec)
#Jop=(q*t/((Np-1)*hbar))*1i*((diag(ones(Np-1,1),-1))-(diag(ones(Np-1,1),1)));

#energy grid
NE = 501
E = np.linspace(- 0.5,1,NE)
dE = E(2) - E(1)
zplus = 1j * 1e-12

pt1 = Ns + 1
pt2 = Ns + Nc
Ho = T(np.array([np.arange(pt1,pt2+1)]),np.array([np.arange(pt1,pt2+1)]))
s_Ho,__ = Ho.shape
tol_outer = 0.0005
tol_inner = 0.001
#VV is the external bias vector
#VV=[0;0.001;0.002;0.003;0.004;0.005;0.006;0.0065;0.007];  #use this to run
VV = 0.2
#multiple voltage bias simulations sequentially
#VV=0;   #Use this for single voltage bias values
#VV=0.004:0.0005:0.0075;
VV = np.transpose(VV)

#defined.
__,sz_V = VV.shape
sum_qc = np.zeros((1,sz_V))
II = np.zeros((1,sz_V))
## Simulation cell. Refer to thesis for details on the algorithm.
for iV in np.arange(1,sz_V+1).reshape(-1):
    flag = 0
    flag2 = 0
    U_guess = np.array([[0 * np.ones((Ns,1))],[0 * np.ones((Nc,1))],[0 * np.ones((Ns,1))]])
    Ec_tot = Ec + U_guess
    conv_outer = 1
    V = VV(iV)
    mu1 = mu - V
    mu2 = mu
    f1 = n0H * np.log(1 + np.exp((mu1 - E) / kTH))
    f2 = n0C * np.log(1 + np.exp((mu2 - E) / kTC))
    for sz in np.arange(1,pt1 - 1+1).reshape(-1):
        N1[sz] = N_cH * Fhalf(- (Ec_tot(sz) - mu1) / kTH)
    sig1 = np.zeros((s_Ho,s_Ho))
    sig2 = np.zeros((s_Ho,s_Ho))
    rho = np.zeros((s_Ho,s_Ho))
    for k in np.arange(1,NE+1).reshape(-1):
        ck = 1 - ((E(k) + zplus - Ec_tot(pt1)) / (2 * t))
        ka = np.arccos(ck)
        sig1[1,1] = - t * np.exp(1j * ka)
        gam1 = 1j * (sig1 - np.transpose(sig1))
        ck = 1 - ((E(k) + zplus - Ec_tot(pt2)) / (2 * t))
        ka = np.arccos(ck)
        sig2[s_Ho,s_Ho] = - t * np.exp(1j * ka)
        gam2 = 1j * (sig2 - np.transpose(sig2))
        G = inv(((E(k) + zplus) * np.eye(s_Ho)) - Ho - diag(Ec_tot(np.arange(pt1,pt2+1))) - sig1 - sig2)
        A1 = np.transpose(G) * gam1 * G
        A2 = np.transpose(G) * gam2 * G
        rho = rho + (dE * ((f1(k) * A1) + (f2(k) * A2)) / (2 * np.pi))
    N2 = (1 / a) * real(diag(rho))
    ct = 1
    for sz in np.arange(pt2 + 1,Np+1).reshape(-1):
        N3[ct] = N_cC * Fhalf(- (Ec_tot(sz) - mu2) / kTC)
        ct = ct + 1
    N = np.array([[np.transpose(N1)],[N2],[np.transpose(N3)]])
    D = np.zeros((Np,1))
    #Newton-Raphson solution of the Poission equation starts here
    #conv_outer checks for convergence of the Poission potential. If the
#difference between the updated and old solution is below tol_outer, the
#Poission solution is converged.
    #conv_innter checks for convergence of the Newton-Raphson solution for the
#Poission potential, since NR itself is an iterative procedure. So, for
#each iteration of the outer loop, the NR loop iterates until
#conv_inner<tol_inner.
    #Fine gridding starts here. The idea here is to first work with a course energy grid till conv_outer reduces to a
#threshold (larger than tol_outer,1e-2 here), and then refine the grid near the
#transmission peaks.
    while (np.abs(conv_outer) > tol_outer):

        U_old = U_guess
        conv_inner = 1
        if (np.abs(conv_outer) > 0.01 and flag == 0):
            NE = 1001
            E = np.linspace(- 0.5,1,NE)
            dE = E(2) - E(1)
            dE_old = dE
            f1 = n0H * np.log(1 + np.exp((mu1 - E) / kTH))
            f2 = n0C * np.log(1 + np.exp((mu2 - E) / kTC))
        else:
            pks,posn = findpeaks(T12)
            __,sz_posn = posn.shape
            ctr = 0
            for i_posn in np.arange(1,sz_posn+1).reshape(-1):
                if (pks(i_posn) > 1e-05 and ctr == 0):
                    res1 = E(posn(i_posn))
                    ctr = 1
            #res1=E(posn(1));
#res2=E(posn(2));
            Emin = - 0.5
            Emax = 1
            dE_course = dE_old
            dE_fine = 5e-06
            E = np.array([[np.transpose((np.arange(Emin,res1 - 10 * dE_course+dE_course,dE_course)))],[np.transpose((np.arange(res1 - 10 * dE_course,res1 + 10 * dE_course+dE_fine,dE_fine)))],[np.transpose((np.arange(res1 + 10 * dE_course,Emax+dE_course,dE_course)))]])
            NE,__ = E.shape
            #NE=1001;E=linspace(-.25,0.8,NE);dE=E(2)-E(1);
            f1 = n0H * np.log(1 + np.exp((mu1 - E) / kTH))
            f2 = n0C * np.log(1 + np.exp((mu2 - E) / kTC))
            flag = 1
            flag2 = 1
        #Fine gridding ends here
        while (np.abs(conv_inner) > tol_inner):

            x = N / N_c
            Ef = np.log(x) + 0.353553 * x - 0.00495009 * x ** 2 + 0.000148386 * x ** 3 - 4.42563e-06 * x ** 4
            for k in np.arange(1,Np+1).reshape(-1):
                if (Ef(k) < - 10):
                    Ef[k] = - 10
                z = (- Ec_tot(k) + Ef(k) * kT) / kT
                D[k,1] = 2 * ((n0 / 2) ** 1.5) * ((Fhalf(z + 0.001) - Fhalf(z)) / 0.001) / kT
            dN = beta * (N - Nd) + ((1) * D2 * U_guess)
            dU = - (inv(D2 - (beta * diag(D)))) * dN
            U_guess = U_guess + dU
            Ec_tot = Ec + U_guess
            for jk in np.arange(1,Np+1).reshape(-1):
                N[jk] = N_c * Fhalf((- Ec_tot(jk) + Ef(jk) * kT) / kT)
            conv_inner = (np.amax(np.abs(dN))) / (np.amax(np.amax(Nd)))

        for sz in np.arange(1,pt1 - 1+1).reshape(-1):
            N1[sz] = N_cH * Fhalf(- (Ec_tot(sz) - mu1) / kTH)
        sig1 = np.zeros((s_Ho,s_Ho))
        sig2 = np.zeros((s_Ho,s_Ho))
        rho = np.zeros((s_Ho,s_Ho))
        I = 0
        I_h1 = 0
        I_h2 = 0
        for k in np.arange(1,NE+1).reshape(-1):
            ck = 1 - ((E(k) + zplus - Ec_tot(pt1)) / (2 * t))
            ka = np.arccos(ck)
            sig1[1,1] = - t * np.exp(1j * ka)
            gam1 = 1j * (sig1 - np.transpose(sig1))
            ck = 1 - ((E(k) + zplus - Ec_tot(pt2)) / (2 * t))
            ka = np.arccos(ck)
            sig2[s_Ho,s_Ho] = - t * np.exp(1j * ka)
            gam2 = 1j * (sig2 - np.transpose(sig2))
            G = inv(((E(k) + zplus) * np.eye(s_Ho)) - Ho - diag(Ec_tot(np.arange(pt1,pt2+1))) - sig1 - sig2)
            A1 = np.transpose(G) * gam1 * G
            A2 = np.transpose(G) * gam2 * G
            T12[k] = real(trace(gam1 * G * gam2 * np.transpose(G)))
            xH = (mu1 - E(k)) / kTH
            xC = (mu2 - E(k)) / kTC
            tstep = 0.001
            tt = np.arange(0,25+tstep,tstep)
            K2H = n0H * kTH * sum(np.multiply(tt / (1 + np.exp(tt - xH)),tstep))
            K2C = n0C * kTC * sum(np.multiply(tt / (1 + np.exp(tt - xC)),tstep))
            #For k=1 the grid size will definitely be dE_old. for larger k
#however grid size can be dE_old or dE_new (if k is near a peak).
            if k == 1:
                rho = rho + (dE_old * ((f1(k) * A1) + (f2(k) * A2)) / (2 * np.pi))
                I = I + (dE_old * IE * T12(k) * (f1(k) - f2(k)))
                I_h1 = I_h1 + (dE_old * IE * (E(k) - mu1) * T12(k) * (f1(k) - f2(k)))
                I_h2 = I_h2 + (dE_old * IE * T12(k) * (K2H - K2C))
                #I_h2=I_h2+(dE_old*IE*T12(k)*(kTH*f1(k)-kTC*f2(k)));
                I2[k] = IE * T12(k) * (f1(k) - f2(k))
                Ih1[k] = (E(k) - mu1) * IE * T12(k) * (f1(k) - f2(k))
                Ih2[k] = IE * T12(k) * (K2H - K2C)
            else:
                rho = rho + ((E(k) - E(k - 1)) * ((f1(k) * A1) + (f2(k) * A2)) / (2 * np.pi))
                I = I + ((E(k) - E(k - 1)) * IE * T12(k) * (f1(k) - f2(k)))
                I_h1 = I_h1 + ((E(k) - E(k - 1)) * (E(k) - mu1) * IE * T12(k) * (f1(k) - f2(k)))
                I_h2 = I_h2 + ((E(k) - E(k - 1)) * IE * T12(k) * (K2H - K2C))
                #I_h2=I_h2+((E(k)-E(k-1))*IE*T12(k)*(kTH*f1(k)-kTC*f2(k)));
                I2[k] = IE * T12(k) * (f1(k) - f2(k))
                Ih1[k] = (E(k) - mu1) * IE * T12(k) * (f1(k) - f2(k))
                Ih2[k] = IE * T12(k) * (K2H - K2C)
                EE2[k] = E(k) - E(k - 1)
            k
        N2 = (1 / a) * real(diag(rho))
        ct = 1
        for sz in np.arange(pt2 + 1,Np+1).reshape(-1):
            N3[ct] = N_cC * Fhalf(- (Ec_tot(sz) - mu2) / kTC)
            ct = ct + 1
        N = np.array([[np.transpose(N1)],[N2],[np.transpose(N3)]])
        conv_outer = np.amax(np.abs(U_guess - U_old))
        #figure();
#plot(E,T12);
#hold on;
#if (conv_outer<=1e-3)
#    plot(E,T12)
#   hold on
#end

    ##
#semilogy(XX,N);
#hold on;
#plot(XX,Ec_tot);
#hold off;
#N_well=N(Ns+Nl/2+Nb+1:Ns+Nl/2+Nb+Ng);
#sum_qc(iV)=sum(N_well);
    II[iV] = I
    IH1[iV] = I_h1
    IH2[iV] = I_h2
    IH[iV] = IH1(iV) + IH2(iV)
    nu[iV] = II(iV) * VV(iV) / IH(iV)
    #TT(iV,:)=T12;
#EC_TOT(iV,:)=Ec_tot;
#figure();
#semilogy(E,T12);
#hold on;

##
#plot(VV,sum_qc);
#hold on;
##

# I=0;
# sig1=zeros(Np,Np);sig2=zeros(Np,Np);
# rho=0;
# NE=15001;E=linspace(-.25,1,NE);dE=E(2)-E(1);
# f1=n0*log(1+exp((mu-E)./kT));f2=n0*log(1+exp((mu-V-E)./kT));
# for k=1:NE
#     ck=1-((E(k)+zplus-Ec_tot(1))/(2*t));ka=acos(ck);
#     sig1(1,1)=-t*exp(1i*ka);gam1=1i*(sig1-sig1');
#     ck=1-((E(k)+zplus-Ec_tot(Np))/(2*t));ka=acos(ck);
#     sig2(Np,Np)=-t*exp(1i*ka);gam2=1i*(sig2-sig2');
#     G=inv(((E(k)+zplus)*eye(Np))-T-diag(Ec_tot)-sig1-sig2);
#     T12(k)=real(trace(gam1*G'*gam2*G));
#     I=I+(dE*IE*T12(k)*(f1(k)-f2(k)));
# end
##


# ##
#Hxx=T+diag(U_guess);
#[EGN,FN]=eig(Hxx);
#F=sort(F);
#   hold on;
#  plot(VV,nu*11);
toc

def Fhalf(x = None): 
xx = np.linspace(0,np.abs(x) + 10,251)
dx = xx(2) - xx(1)
fx = (2 * dx / np.sqrt(np.pi)) * np.sqrt(xx) / (1 + np.exp(xx - x))
y = sum(fx)
return y
