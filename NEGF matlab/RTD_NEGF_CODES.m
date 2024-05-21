clear;clc;


%% Constant definition cell
hbar=1.06e-34;q=1.6e-19;epsil=10*8.85E-12;
kTH=8.314/6.023e23*330/q;     %kT for hot contact
kTC=8.314/6.023e23*300/q;     %kT for cold contact
kT=0.5*(kTH+kTC);             %Average kT, used later in calculations
m=.07*9.1e-31;n0=2*m*kT*q/(2*pi*(hbar^2));
IE=(q*q)/(2*pi*hbar);
n0H=2*m*kTH*q/(2*pi*(hbar^2));
n0C=2*m*kTC*q/(2*pi*(hbar^2));
tic;
a=2.5e-10;   %Spatial grid step size(a<2.5e-10 can be used, but is < GaAs lattice constant)
t=(hbar^2)/(2*m*(a^2)*q);
beta=q*a*a/epsil;

%% Device definition
%Device is divided into 3 parts: Heavily doped contacts (source and drain),
% lightly doped active region and lightly doped spacer (between contacts
% and active region). The active region consists of the two barriers and
% well. The active region together with the spacers is referred to as the
% channel.

%The "quantum" simulation is only performed in the channel region, since
%confinement only exists here. In the contact regions a semiclassical
%equation is used for electron population, as described later. Quantum
%simulation for the whole device could be used, but it would lead to a much
%larger computational time.

%In principle the effective mass and permittivity should be different for
%barrier (AlAs) and well (GaAs). however since the electron concentration
%in the barrier is really small, and the permittivities of GaAs and AlAs
%are similar, both effective mass and permittivity are taken to be constant
%throughout the device. This should, however, be corrected in the long
%run, and would require manual summing over all lateral modes (with
%constant effective mass throughout, this summing results in the F2D
%functions being used for electron population instead of the fermi
%function). I have done this for a simpler problem (self-consistent
%Schrodinger-Poission solver); please contact me in case manual sum over
%modes proves problematic.
%Only a small effect is expected on calculations.
Ls=30e-9;                   %Length of source/drain  
Lbl=4e-9;                   %length of left barrier
Lbr=4e-9;                   %length of right barrier
Lg=7e-9;                    %length of well
Lbuff=1e-9;                 %length of buffer/spacer region (between heavily doped source and lightly doped active region)
Ns=round(Ls/a);
Nbl=round(Lbl/a);
Nbr=round(Lbr/a);
Ng=round(Lg/a);
Nl=round(2*Lbuff/a);  %Buffer is present between source and active region as well as drain and active region
Nc=Nl+Nbl+Nbr+Ng;    %Total size of "channel" (Channel=active region + spacer)
Np=Ns+Nc+Ns;XX=a*1e9*[1:1:Np];
mu=0;
N_cH=2*((n0H/2)^1.5); %Effective density of states near hot contact
N_cC=2*((n0C/2)^1.5); %effective density of states near cold contact
N_c=2*((n0/2)^1.5);  %Average effective density of states
Nd=1e24;   %Units of /m^3
Nd=Nd*[ones(Ns,1);.001*ones(Nc,1);ones(Ns,1)];  %Doping profile. Heavily doped contacts, lightly doped channel

D2=-(2*diag(ones(1,Np)))+(diag(ones(1,Np-1),1))+(diag(ones(1,Np-1),-1));  %Tridiagonal second derivative matrix
D2(1,1)=-1;D2(Np,Np)=-1;%zero field (Neumann) boundary condition enforced by setting (1,1) and (Np,Np) elements to -1 instead of -2 (-2 would be for Dirichelt)
T=(2*t*diag(ones(1,Np)))-(t*diag(ones(1,Np-1),1))-(t*diag(ones(1,Np-1),-1));
Ec=[zeros(Ns,1);zeros(Nl/2,1);0.1*ones(Nbl,1);zeros(Ng,1);0.1*ones(Nbr,1);zeros(Nl/2,1);zeros(Ns,1)];
T=T+diag(Ec);
%Jop=(q*t/((Np-1)*hbar))*1i*((diag(ones(Np-1,1),-1))-(diag(ones(Np-1,1),1)));

%energy grid
NE=501;E=linspace(-0.5,1,NE);dE=E(2)-E(1);zplus=1i*1e-12;  %NE is the energy grid size, dE is the step

pt1=Ns+1;
pt2=Ns+Nc;
Ho=T([pt1:pt2],[pt1:pt2]);
[s_Ho,~]=size(Ho);
tol_outer=5e-4;
tol_inner=1e-3;

%VV is the external bias vector
%VV=[0;0.001;0.002;0.003;0.004;0.005;0.006;0.0065;0.007];  %use this to run
VV=0.2;
%multiple voltage bias simulations sequentially
%VV=0;   %Use this for single voltage bias values
%VV=0.004:0.0005:0.0075;
VV=VV';%This may or may not be required, based on how the VV vector is
%defined.
[~,sz_V]=size(VV);
sum_qc=zeros(1,sz_V);
II=zeros(1,sz_V);
%% Simulation cell. Refer to thesis for details on the algorithm. 
for iV=1:sz_V
    flag=0;
    flag2=0;
    U_guess=[0*ones(Ns,1);0*ones(Nc,1);0*ones(Ns,1)];  %initial guess for the poission potential
    Ec_tot=Ec+U_guess;  %Add the conduction band profile to U_guess to get the final profile of the CB
    conv_outer=1;
    V=VV(iV)
    mu1=mu-V;
    mu2=mu;
    f1=n0H*log(1+exp((mu1-E)./kTH));f2=n0C*log(1+exp((mu2-E)./kTC));  %f1 and f2 are the f2D functions at the two ends of the device
for sz=1:pt1-1
    N1(sz)=N_cH*Fhalf(-(Ec_tot(sz)-mu1)/kTH);   %Semiclassical expression for electron concentration in the source region 
end
sig1=zeros(s_Ho);sig2=zeros(s_Ho);rho=zeros(s_Ho);
for k=1:NE
    ck=1-((E(k)+zplus-Ec_tot(pt1))/(2*t));ka=acos(ck);
    sig1(1,1)=-t*exp(1i*ka);gam1=1i*(sig1-sig1');
    ck=1-((E(k)+zplus-Ec_tot(pt2))/(2*t));ka=acos(ck);
    sig2(s_Ho,s_Ho)=-t*exp(1i*ka);gam2=1i*(sig2-sig2');
    G=inv(((E(k)+zplus)*eye(s_Ho))-Ho-diag(Ec_tot(pt1:pt2))-sig1-sig2);
    A1=G'*gam1*G;A2=G'*gam2*G;  
    rho=rho+(dE*((f1(k)*A1)+(f2(k)*A2))/(2*pi));
end
N2=(1/a)*real(diag(rho));     %NEGF expression for electron density in channel
ct=1;
for sz=pt2+1:Np
    N3(ct)=N_cC*Fhalf(-(Ec_tot(sz)-mu2)/kTC);  %Semiclassical expression for electron concentration in the source region 
    ct=ct+1;
end
N=[N1';N2;N3'];
D=zeros(Np,1);
%Newton-Raphson solution of the Poission equation starts here

%conv_outer checks for convergence of the Poission potential. If the
%difference between the updated and old solution is below tol_outer, the
%Poission solution is converged.

%conv_innter checks for convergence of the Newton-Raphson solution for the
%Poission potential, since NR itself is an iterative procedure. So, for
%each iteration of the outer loop, the NR loop iterates until
%conv_inner<tol_inner.

%Fine gridding starts here. The idea here is to first work with a course energy grid till conv_outer reduces to a
%threshold (larger than tol_outer,1e-2 here), and then refine the grid near the
%transmission peaks. 

while(abs(conv_outer)>tol_outer)
    U_old=U_guess;
    conv_inner=1;
    if (abs(conv_outer)>1e-2 && flag==0)
        NE=1001;E=linspace(-0.5,1,NE);dE=E(2)-E(1);  %Initial Course grid 
        dE_old=dE;
        f1=n0H*log(1+exp((mu1-E)./kTH));f2=n0C*log(1+exp((mu2-E)./kTC));    
    else 
        [pks,posn]=findpeaks(T12);  %Finding transmission peaks
        [~,sz_posn]=size(posn);
        ctr=0;
        for i_posn=1:sz_posn
            if(pks(i_posn)>1e-5&&ctr==0)   %Checking if peaks are significant. Sometimes there can be spurious spikes in the transmission function which need to be ignored, hence a threshold of 1e-5 is set
                res1=E(posn(i_posn));
                ctr=1;
            end
        end
        %res1=E(posn(1));
        %res2=E(posn(2));
        Emin=-0.5;
        Emax=1;
        dE_course=dE_old; %The grid size is same as the old grid, except near the ground state peak where dE_fine is used
        dE_fine=5e-6;  %Fine grid size
        E=[(Emin:dE_course:res1-10*dE_course)';(res1-10*dE_course:dE_fine:res1+10*dE_course)';(res1+10*dE_course:dE_course:Emax)'];  %Fine grid around first transmission peak only, since that is the peak which is difficult to resolve. The fine grid is applied till 10*dE_course on either side of the peak. 10 is arbitrary, can play around to optimize.
        [NE,~]=size(E);
        %NE=1001;E=linspace(-.25,0.8,NE);dE=E(2)-E(1);
        f1=n0H*log(1+exp((mu1-E)./kTH));f2=n0C*log(1+exp((mu2-E)./kTC));
        flag=1;
        flag2=1;
    end
%Fine gridding ends here

    while(abs(conv_inner)>tol_inner)
        x=N/N_c;
        Ef = log(x) + 3.53553e-1*x - 4.95009e-3*x.^2 + 1.48386e-4*x.^3 - 4.42563e-6*x.^4;   %Extracting a pseudo fermi level from the electron population. This is to be used in the NR solution
        for k=1:Np
            if(Ef(k)<-10)
                Ef(k)=-10;    %Ef can shoot off to really low values near the barrier (since the leading term is log (x)), which messes up calculations. Thus a lowe bound is put on it.
            end
            z=(-Ec_tot(k)+Ef(k)*kT)/kT;  
            D(k,1)=2*((n0/2)^1.5)*((Fhalf(z+.001)-Fhalf(z))/.001)/kT; 
        end
        dN=beta*(N-Nd)+((1)*D2*U_guess);     %Refer the appendix of Lake and Klimeck's paper on RTD
        dU=-(inv(D2-(beta*diag(D))))*dN;U_guess=U_guess+dU;
        Ec_tot=Ec+U_guess;
        for jk=1:Np
            N(jk)=N_c*Fhalf((-Ec_tot(jk)+Ef(jk)*kT)/kT);
        end
        conv_inner=(max(abs(dN)))/(max(max(Nd)));
    end
    for sz=1:pt1-1
        N1(sz)=N_cH*Fhalf(-(Ec_tot(sz)-mu1)/kTH);
    end
    sig1=zeros(s_Ho);sig2=zeros(s_Ho);rho=zeros(s_Ho);
    I=0;
    I_h1=0;
    I_h2=0;
    for k=1:NE
        ck=1-((E(k)+zplus-Ec_tot(pt1))/(2*t));ka=acos(ck);
        sig1(1,1)=-t*exp(1i*ka);gam1=1i*(sig1-sig1');
        ck=1-((E(k)+zplus-Ec_tot(pt2))/(2*t));ka=acos(ck);
        sig2(s_Ho,s_Ho)=-t*exp(1i*ka);gam2=1i*(sig2-sig2');
        G=inv(((E(k)+zplus)*eye(s_Ho))-Ho-diag(Ec_tot(pt1:pt2))-sig1-sig2);
        A1=G'*gam1*G;A2=G'*gam2*G;  
        T12(k)=real(trace(gam1*G*gam2*G'));
        xH=(mu1-E(k))/kTH;
        xC=(mu2-E(k))/kTC;
        tstep=0.001;
        tt=0:tstep:25;
        K2H=n0H*kTH*sum(tt./(1+exp(tt-xH)).*tstep);   %Refer my APL paper on RTD based TE for clarification on K2, Ih1 and Ih2.
        K2C=n0C*kTC*sum(tt./(1+exp(tt-xC)).*tstep);
        %For k=1 the grid size will definitely be dE_old. for larger k
        %however grid size can be dE_old or dE_new (if k is near a peak).
        if k==1;
            rho=rho+(dE_old*((f1(k)*A1)+(f2(k)*A2))/(2*pi));
            I=I+(dE_old*IE*T12(k)*(f1(k)-f2(k)));
            I_h1=I_h1+(dE_old*IE*(E(k)-mu1)*T12(k)*(f1(k)-f2(k)));
            I_h2=I_h2+(dE_old*IE*T12(k)*(K2H-K2C));
            %I_h2=I_h2+(dE_old*IE*T12(k)*(kTH*f1(k)-kTC*f2(k)));
            I2(k)=IE*T12(k)*(f1(k)-f2(k));
            Ih1(k)=(E(k)-mu1)*IE*T12(k)*(f1(k)-f2(k));
            Ih2(k)=IE*T12(k)*(K2H-K2C);
            
        else
            rho=rho+((E(k)-E(k-1))*((f1(k)*A1)+(f2(k)*A2))/(2*pi));   %E(k)-E(k-1) takes care of varying grid sizes.
            I=I+((E(k)-E(k-1))*IE*T12(k)*(f1(k)-f2(k)));
            I_h1=I_h1+((E(k)-E(k-1))*(E(k)-mu1)*IE*T12(k)*(f1(k)-f2(k)));
            I_h2=I_h2+((E(k)-E(k-1))*IE*T12(k)*(K2H-K2C));
            %I_h2=I_h2+((E(k)-E(k-1))*IE*T12(k)*(kTH*f1(k)-kTC*f2(k)));
            I2(k)=IE*T12(k)*(f1(k)-f2(k));
            Ih1(k)=(E(k)-mu1)*IE*T12(k)*(f1(k)-f2(k));
            Ih2(k)=IE*T12(k)*(K2H-K2C);
            EE2(k)=E(k)-E(k-1);
        end
        k;
    end
    N2=(1/a)*real(diag(rho));
    ct=1;
    for sz=pt2+1:Np
        N3(ct)=N_cC*Fhalf(-(Ec_tot(sz)-mu2)/kTC);
        ct=ct+1;
    end
    N=[N1';N2;N3'];
    conv_outer=max(abs(U_guess-U_old));
    %figure();
    %plot(E,T12);
    %hold on;
    %if (conv_outer<=1e-3)
    %    plot(E,T12)
    %   hold on
    %end
end
%%  
%semilogy(XX,N);
%hold on;
%plot(XX,Ec_tot);
%hold off;
%N_well=N(Ns+Nl/2+Nb+1:Ns+Nl/2+Nb+Ng);
%sum_qc(iV)=sum(N_well);
II(iV)=I;
IH1(iV)=I_h1;
IH2(iV)=I_h2;
IH(iV)=IH1(iV)+IH2(iV);
nu(iV)=II(iV)*VV(iV)/IH(iV);
%TT(iV,:)=T12;
%EC_TOT(iV,:)=Ec_tot;
%figure();
%semilogy(E,T12);
%hold on;
end
%%
%plot(VV,sum_qc);
%hold on;
%%
% 
% I=0;
% sig1=zeros(Np,Np);sig2=zeros(Np,Np);
% rho=0;
% NE=15001;E=linspace(-.25,1,NE);dE=E(2)-E(1);
% f1=n0*log(1+exp((mu-E)./kT));f2=n0*log(1+exp((mu-V-E)./kT));
% for k=1:NE
%     ck=1-((E(k)+zplus-Ec_tot(1))/(2*t));ka=acos(ck);
%     sig1(1,1)=-t*exp(1i*ka);gam1=1i*(sig1-sig1');
%     ck=1-((E(k)+zplus-Ec_tot(Np))/(2*t));ka=acos(ck);
%     sig2(Np,Np)=-t*exp(1i*ka);gam2=1i*(sig2-sig2');
%     G=inv(((E(k)+zplus)*eye(Np))-T-diag(Ec_tot)-sig1-sig2);  
%     T12(k)=real(trace(gam1*G'*gam2*G));
%     I=I+(dE*IE*T12(k)*(f1(k)-f2(k)));
% end
%%

% 
% %%
 %Hxx=T+diag(U_guess);
 %[EGN,FN]=eig(Hxx);
 %F=sort(F);
%   hold on;
%  plot(VV,nu*11);
toc;

function y=Fhalf(x)
xx=linspace(0,abs(x)+10,251);dx=xx(2)-xx(1);
fx=(2*dx/sqrt(pi))*sqrt(xx)./(1+exp(xx-x));
y=sum(fx);
end
