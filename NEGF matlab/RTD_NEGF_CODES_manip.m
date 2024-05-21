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

%%%%%Time start%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Change length as needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ls=2.45e-8;                   %Length of source/drain  
Lbl=1.2e-8;                   %length of left barrier
Lbr=1.2e-8;                   %length of right barrier
Lg=6e-9;                      %length of well
Lbuff=5e-10;              %length of buffer/spacer region (between heavily doped source and lightly doped active region)
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
Nd=1e25;   %Units of /m^3
Nd=Nd*[ones(Ns,1);0.001*ones(Nc,1);ones(Ns,1)];  %Doping profile. Heavily doped contacts, lightly doped channel


D2=-(2*diag(ones(1,Np)))+(diag(ones(1,Np-1),1))+(diag(ones(1,Np-1),-1));  %Tridiagonal second derivative matrix
D2(1,1)=-1;D2(Np,Np)=-1;            %zero field (Neumann) boundary condition enforced by setting (1,1) and (Np,Np) elements to -1 instead of -2 (-2 would be for Dirichelt)
T=(2*t*diag(ones(1,Np)))-(t*diag(ones(1,Np-1),1))-(t*diag(ones(1,Np-1),-1));
Ec=[zeros(Ns,1);zeros(Nl/2,1);0.1*ones(Nbl,1);zeros(Ng,1);0.1*ones(Nbr,1);zeros(Nl/2,1);zeros(Ns,1)];
T=T+diag(Ec);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Change bias as needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VV=0.0;

VV=VV';%This may or may not be required, based on how the VV vector is
%defined.
[~,sz_V]=size(VV);
sum_qc=zeros(1,sz_V);
II=zeros(1,sz_V);
%% Simulation cell. Refer to thesis for details on the algorithm. 
for iV=1:sz_V
    flag=0;
    flag2=0;
    %%%%%%%%%% As per the PP - FB value from model %%%%%%%%%%%%%
    %U_guess=[0*ones(Ns,1);0*ones(Nc,1);0*ones(Ns,1)];  %initial guess for the poission potential
    U_guess=[-0.24187005;-0.24191728;-0.24192815;-0.24187462;-0.24191032;-0.24189034;-0.24202852;-0.24199781;-0.24199797;-0.24190862;-0.24190865;-0.24184152;-0.24197346;-0.24181716;-0.241956;-0.24188276;-0.24197403;-0.24181515;-0.24181563;-0.24187498;-0.24208875;-0.24183223;-0.24175592;-0.24179055;-0.24186447;-0.24184088;-0.24185686;-0.24191013;-0.24192064;-0.24179545;-0.24179593;-0.24184023;-0.24178325;-0.24193266;-0.2420169;-0.24189176;-0.2417707;-0.24188389;-0.24197987;-0.24184655;-0.24180108;-0.24192227;-0.24177444;-0.2419326;-0.24170238;-0.24180886;-0.24189268;-0.24186933;-0.24185705;-0.24181885;-0.24181183;-0.24187146;-0.24046366;-0.24173501;-0.24080628;-0.24060936;-0.24047604;-0.24189554;-0.24033928;-0.2417507;-0.24078417;-0.24052358;-0.24192007;-0.24072343;-0.24085343;-0.2407238;-0.24067175;-0.24047968;-0.2410433;-0.24036503;-0.24112667;-0.24488315;-0.24456984;-0.24417861;-0.24447039;-0.2442279;-0.24455509;-0.24448688;-0.244538;-0.24459797;-0.24348159;-0.24325737;-0.24267104;-0.24256757;-0.24207853;-0.24343507;-0.24609949;-0.24619865;-0.24785362;-0.24605367;-0.24789701;-0.2453358;-0.24724893;-0.24481504;-0.24779217;-0.24870618;-0.2349497;-0.23939456;-0.23931263;-0.24156524;-0.3283241;-0.32971582;-0.32150713;-0.31300455;-0.31202275;-0.30104953;-0.30267024;-0.29375392;-0.28986576;-0.2887406;-0.22549629;-0.23521166;-0.2137167;-0.21873443;-0.21222107;-0.20637682;-0.1414586;-0.14357601;-0.16071367;-0.15382382;-0.13763462;-0.13676438;-0.17267047;-0.17791763;-0.13058068;-0.14058919;-0.1712511;-0.16687064;-0.20358594;-0.18979602;-0.26348737;-0.27128604;-0.27585632;-0.2618736;-0.2863672;-0.2757555;-0.33215204;-0.3215759;-0.30178896;-0.31346023;-0.30903837;-0.30500767;-0.30412444;-0.30195537;-0.29890501;-0.29484943;-0.21564187;-0.18142757;-0.1412921;-0.1547929;-0.14257526;-0.16015281;-0.10667794;-0.07893988;-0.09766509;-0.097161;-0.10654436;-0.105675675;-0.14914533;-0.1430553;-0.09705501;-0.09148721;-0.1848485;-0.19888283;-0.20249467;-0.17089461;-0.18614355;-0.17172407;-0.15578094;-0.14185525;-0.12896736;-0.13381787;-0.28917968;-0.27987105;-0.28156275;-0.28873423;-0.2909766;-0.28670126;-0.21650267;-0.23018233;-0.22575185;-0.24416931;-0.20501263;-0.19491954;-0.14432032;-0.159899;-0.12194187;-0.10273867;-0.13725731;-0.15310277;-0.14998141;-0.14489257;-0.07961282;-0.1082533;-0.19173968;-0.19354495;-0.18739295;-0.18422407;-0.19938299;-0.19485256;-0.20227386;-0.20596653;-0.22984582;-0.23164007;-0.23460266;-0.23926884;-0.26644227;-0.26119596;-0.26477548;-0.26836538;-0.21807484;-0.20302962;-0.23371239;-0.19574031;-0.19344372;-0.18079287;-0.22918978;-0.2207995;-0.17609066;-0.17059195;-0.106265664;-0.13159084;-0.11042634;-0.10184595;-0.075936906;-0.05847508;-0.10355552;-0.095840134;-0.12527168;-0.12307425;-0.20765539;-0.20472446;-0.2161362;-0.21523546;-0.2001978;-0.20540203;-0.20897512;-0.20732261;-0.20602529;-0.20670609;-0.23355396;-0.21758042;-0.24945003;-0.24566579;-0.24282272;-0.24301226;-0.23936346;-0.24487776;-0.24039888;-0.24068774;-0.24111938;-0.24281693;-0.23765503;-0.23779441;-0.25066745;-0.23309956;-0.23781952;-0.23335063;-0.23310839;-0.2337055;-0.23349248;-0.23351014;-0.22915308;-0.22903061;-0.22914523;-0.22903702;-0.24178982;-0.24244015;-0.2421796;-0.24238743;-0.24232861;-0.24219273;-0.24243595;-0.24245973;-0.24191613;-0.24178614;-0.24232948;-0.24217847;-0.24226956;-0.24247092;-0.2428222;-0.24288645;-0.24294645;-0.24291258;-0.24291691;-0.24291147;-0.24289526;-0.24291775;-0.24291421;-0.24294479;-0.24294741;-0.24251167;-0.24259815;-0.24230854;-0.2424783;-0.24247995;-0.24248144;-0.2424827;-0.24237558;-0.2424852;-0.24245726;-0.24251209;-0.24244942;-0.24245971;-0.24254915;-0.24246137;-0.24251616;-0.24246272;-0.24262457;-0.2425504;-0.24254218;-0.24262603;-0.24239476;-0.2424051;-0.24270694;-0.24254405;-0.24258684;-0.24257383;-0.24257408;-0.24257404];  %initial guess for the poission potential
    Ec_tot=Ec+U_guess;                                 %Add the conduction band profile to U_guess to get the final profile of the CB
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

        Emin=-0.5;
        Emax=1;
        dE_course=dE_old; %The grid size is same as the old grid, except near the ground state peak where dE_fine is used
        dE_fine=5e-6;     %Fine grid size
        E=[(Emin:dE_course:res1-10*dE_course)';(res1-10*dE_course:dE_fine:res1+10*dE_course)';(res1+10*dE_course:dE_course:Emax)'];  %Fine grid around first transmission peak only, since that is the peak which is difficult to resolve. The fine grid is applied till 10*dE_course on either side of the peak. 10 is arbitrary, can play around to optimize.
        [NE,~]=size(E);

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
            I2(k)=IE*T12(k)*(f1(k)-f2(k));
            Ih1(k)=(E(k)-mu1)*IE*T12(k)*(f1(k)-f2(k));
            Ih2(k)=IE*T12(k)*(K2H-K2C);
            
        else
            rho=rho+((E(k)-E(k-1))*((f1(k)*A1)+(f2(k)*A2))/(2*pi));   %E(k)-E(k-1) takes care of varying grid sizes.
            I=I+((E(k)-E(k-1))*IE*T12(k)*(f1(k)-f2(k)));
            I_h1=I_h1+((E(k)-E(k-1))*(E(k)-mu1)*IE*T12(k)*(f1(k)-f2(k)));
            I_h2=I_h2+((E(k)-E(k-1))*IE*T12(k)*(K2H-K2C));
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
    %%%%%%Figure stuff%%%%%%%%
    %figure();
    %plot(E,T12);
    %hold on;
    %if (conv_outer<=1e-3)
    %    plot(E,T12)
    %   hold on
    %end

end
%%
%%%%%%Figure stuff%%%%%%%%
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

%%%%%%Figure stuff%%%%%%%%
%TT(iV,:)=T12;
%EC_TOT(iV,:)=Ec_tot;
%figure();
%semilogy(E,T12);
%hold on;
end
%%

%%%%%%Figure stuff%%%%%%%%
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
%%%%%%%Time ends%%%%%%%%

figure();
plot(Ec_tot,'color','r','linewidth',4);

function y=Fhalf(x)
xx=linspace(0,abs(x)+10,251);dx=xx(2)-xx(1);
fx=(2*dx/sqrt(pi))*sqrt(xx)./(1+exp(xx-x));
y=sum(fx);
end
