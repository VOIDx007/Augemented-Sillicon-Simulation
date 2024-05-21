clear;clc;

a = 2.5e-10;

%Calculate the value of Ns from the starting Doping value of doping
%profile
%Calculate the value of Nl/2 from the flat band profile which will be the
%initial trail of zeros minus the value of Ns
%Calculate Nbl and Nbr as the number of 0.1 in the flat band profile
%Get Ng similarly as the contained zeros between the Nbl and Nbr 0.1s
%Check the value of Doping profile in between and change the value of Nd
%matrix accordingly
%Change the bias voltage appropriately
%First non converging lengths%
% Bias = 0.0;
% Nd = 1e25;
% Nd=Nd*[ones(Ns,1);0.001*ones(Nc,1);ones(Ns,1)];

% Ls=2.45e-8;                   %Length of source/drain  
% Lbl=1.2e-8;                   %length of left barrier
% Lbr=1.2e-8;                   %length of right barrier
% Lg=6e-9;                      %length of well
% Lbuff=5e-10;                  %length of buffer/spacer region (between heavily doped source and lightly doped active region)

%Second non converging lengths%
% Bias = 0.0;
% Nd = 1e25;
% Nd=Nd*[ones(Ns,1);0.001*ones(Nc,1);ones(Ns,1)];

% Ls=2.35e-8;                   %Length of source/drain  
% Lbl=1.2e-8;                   %length of left barrier
% Lbr=1.2e-8;                   %length of right barrier
% Lg=8e-9;                      %length of well
% Lbuff=5e-10;                  %length of buffer/spacer region (between heavily doped source and lightly doped active region)

%Third non converging lengths%
% Bias = 0.4;
% Nd = 1e24;
% Nd=Nd*[ones(Ns,1);5*0.001*ones(Nc,1);ones(Ns,1)];

% Ls=2.95e-8;                   %Length of source/drain  
% Lbl=8e-9;                   %length of left barrier
% Lbr=4e-9;                   %length of right barrier
% Lg=8e-9;                      %length of well
% Lbuff=5e-10;

%Fourth non converging lengths%
% Bias = 0.4;
% Nd = 1e24;
% Nd=Nd*[ones(Ns,1);5*0.001*ones(Nc,1);ones(Ns,1)];

% Ls=2.75e-8;                   %Length of source/drain  
% Lbl=8e-9;                   %length of left barrier
% Lbr=8e-9;                   %length of right barrier
% Lg=8e-9;                      %length of well
% Lbuff=5e-10;

%Fifth non converging lengths%
% Bias = 0.3;
% Nd = 1e24;
% Nd=Nd*[ones(Ns,1);0.0001*ones(Nc,1);ones(Ns,1)];

% Ls=2.85e-8;                   %Length of source/drain  
% Lbl=8e-9;                   %length of left barrier
% Lbr=8e-9;                   %length of right barrier
% Lg=6e-9;                      %length of well
% Lbuff=5e-10;

%Sixth non converging lengths%
% Bias = 0.3;
% Nd = 5e24;
% Nd=Nd*[ones(Ns,1);0.001*ones(Nc,1);ones(Ns,1)];

% Ls=2.75e-8;                   %Length of source/drain  
% Lbl=4e-9;                   %length of left barrier
% Lbr=1.2e-8;                   %length of right barrier
% Lg=8e-9;                      %length of well
% Lbuff=5e-10;

%Test Case for values
Ns=round(Ls/a);
Nbl=round(Lbl/a);
Nbr=round(Lbr/a);
Ng=round(Lg/a);
Nl=round(2*Lbuff/a);  %Buffer is present between source and active region as well as drain and active region
Nc=Nl+Nbl+Nbr+Ng;    %Total size of "channel" (Channel=active region + spacer)
Np=Ns+Nc+Ns;XX=a*1e9*[1:1:Np];
mu=0;
fprintf("%d %d %d %d %d %d",Ns, Nbl, Nbr, Ng, Nl, Nc)


Nd=1e25;   %Units of /m^3
Nd=Nd*[ones(Ns,1);5*0.001*ones(Nc,1);ones(Ns,1)];
size(Nd)
Ec=[zeros(Ns,1);zeros(Nl/2,1);0.1*ones(Nbl,1);zeros(Ng,1);0.1*ones(Nbr,1);zeros(Nl/2,1);zeros(Ns,1)];
size(Ec)