% function [X, IncidenceTOT, PARSALL220_79, PARSALL2ge20_79, PAROALL2ge20_79,PAROALL220_79, PARFALL220_79, PARFALL2ge20_79,PrevGeSALL20_79,PrevGeOALL20_79,PrevGeFALL20_79,PrevGeS2079,PrevGeO2079,PrevGeF2079,FALL,OALL,SALL,DMGE2079, PrevDMGE2079, PrevDMALL20_79, S1sumF, S1sumM,PrevDMALL15_64,PrevDMGE,IncidenceRFGE15_64,Nagetot,DMGE,DMNGE,IncidenceRFGE15_64sex,DMAGE2079,IncRate,IncRatesex, PrevOALL, PrevGeO, PrevSALL, PrevGeS, PrevFALL, PrevGeF, PAROALL2, PAROALL2ge, PARSALL2, PARSALL2ge, PARFALL2, PARFALL2ge]=DM_model_3RFs_FINAL(fittingparam1, fittingparam2, data1, BFf, BMm,Btot)
%
%DM MODEL  September 07, 2015
% SUSANNE AWAD
%%
% clear all
% close all
% clc
% load('pram1.mat')
close all
clc 
%% READ XL SHEET
Data1=xlsread('data1-Qatar','DATAALL');
data1=Data1(:,1);   % DEFINE THE DATA MATRIX

PoulationSize='Population Size Zoom.xlsx';
B1ALL=xlsread('Population Size Zoom.xlsx', 'QATAR 1');
B2ALLn=xlsread('Population Size Zoom.xlsx', 'QATAR 2');
B2ALL(1,:,:)=B2ALLn(1:5,:);
B2ALL(2,:,:)=B2ALLn(7:11,:);
BFf=squeeze(B2ALL(1,:,:));%xlsread('Population Size Zoom.xlsx', 'Females');
BMm=squeeze(B2ALL(1,:,:));%xlsread('Population Size Zoom.xlsx', 'Males');
Btot=squeeze(B1ALL(:,:));%xlsread('Population Size Zoom.xlsx', 'sheet3');
%% READ PARAMETER
load('pram1DEC17')       % DM and RF initial info
load('pramOBDEC17')       % DM and RF initial info
load('pramDMDEC17')       % DM and RF initial info
load('ResultQATARnew.mat')  % demog initial info
Result=ResultAll;
%% INITIAL GUESS PARAMETER DEFINITION
% from healthy to Smoker}
beta=param(:,:,1);

% from healthy to PHYSICAL INACTIVE}
shei=param(:,:,2);

%from obese to healthy
sigma=param(:,:,3);

%force of infection
% lambdaH=param(:,:,3);%lambdaH';
%%  FITTING PARAMETER
param_start(:,:,1)=beta;
param_start(:,:,2)=shei;
param_start(:,:,3)=sigma;
fittingparam1=reshape(param_start,20*2*3,1);

% DM incidence rate measures
fittingparam1(121:126,1)=paramDM(:,1);
Estimates1=paramDM(:,1);
% obesity measures
fittingparam1(127:132,1)=paramOB(:,1);
Estimates3=paramOB(:,1);

%demographic measures
fittingparam2=Result;
%%
%%
global ag gen NcompTOT N0 rf NrfTOT delta phi epsilon iota q tt d a rho theta pii ii p ss b cc nu xi sampi w n hh e etaa omega kappa xx ff gg oo gamma upsilon r u k iii Leps chi psi Omg m jj v l 
%%
opengl hardware
set(0,'DefaultLineLineSmoothing','on');
set(0,'DefaultPatchLineSmoothing','on');
opengl('OpenGLLineSmoothingBug',1);
%%
gen=2;
% frac1=1.155;%0.1; % (for Japan) %0.5; $for UK; % to reduce Obesity prevalence so it is equal that of Japan (i.e. 5%)
% frac2=1.455; %physical inactivity
% define input parameter based on the fitting previously done
param1s=fittingparam1(1:20*2*3,1);
param1=reshape(param1s,20,2,3);
beta=param1(:,:,1);
%rates for Smoking and PIA
shei=param1(:,:,2);
sigma=param1(:,:,3);
%DM incidence rate
param3s=fittingparam1(20*2*3+1:(20*2*3)+(6),1,1);
param3=reshape(param3s,3,2,1);
lambdaHh=param3;
lambdaH=(DMincFuncGaussian(ag,gen,lambdaHh(1,:),lambdaHh(2,:),lambdaHh(3,:)));
% figure(1); plot(1:20,lambdaH)
%define rate of beseity
param2s=fittingparam1((20*2*3)+(6)+1:(20*2*3+6)+(3*2),1);
param2=reshape(param2s,3,2,1);
alphas=param2;

% %Population
Result=fittingparam2;

%%
age_gr=5;       % 5-year age stratification of the age groups
ag=100/age_gr;   % Number of age groups
Dsys=100;       % The maximum age included
eta(1:ag,1)=ag/Dsys;       % ageing rate
eta(ag+1,1)=0;
% eta(ag,1)=0;
%%
gnd=2;          % Gender; 1:Femal, 2: uncircumcised men, 3: Circumcised men
rf=3;         % Number of risk factors (from healthy to obese, smoker and physical inactivity)
NrfTOT=2^(rf)-1; % total number of risk factor with overlapping
NcompTOT=NrfTOT*2+2; %total number of compartment that the population is stratified by (including healthy)
%% 
scenario=6.5; % choose between the different scenario I create (1,2,3)(1: Only cohort effect, 3: combination of 2
%                                                             2.sc 4: Increase in Obesity as per WHO (http://gamapserver.who.int/gho/interactive_charts/ncd/risk_factors/obesity/atlas.html), 
%                                                             2. sc 2: increase in Smoking (http://apps.who.int/iris/bitstream/10665/156262/1/9789241564922_eng.pdf )
%                                                             2. sc 3: Increase in Obesity as per WHO (http://gamapserver.who.int/gho/interactive_charts/ncd/risk_factors/obesity/atlas.html), )
intv=2;
IE=0;
% da=zeros(NcompTOT*intv,ag,gen);%zeros(ag,gen);%
% da(1:16,1:20,:)=1;

%%
%% %Time scale
% {use only after first run}
t0=1904;         %Start time
% t0=2010;         %Start time
tf=2050;         %Stop time
% tf=2030;         %Stop time
dt=0.5;         %time interval
tspan=t0:dt:tf;  %timespan to use in the ode45
%% SURVEY DATES
crossdate1=find(tspan==2012); %NATIONAL SURVEY
crossdate2=find(tspan==2008);%NATIONAL SURVEY
crossdate3=find(tspan==2000);%NATIONAL SURVEY
crossdate4=find(tspan==1991);%NATIONAL SURVEY
crossdate5=find(tspan==2010);%SUBNATIONAL SURVEY
crossdate6=find(tspan==2006);%SUBNATIONAL SURVEY
crossdate7=find(tspan==2001);%SUBNATIONAL SURVEY
crossdate8=find(tspan==2030);%SUBNATIONAL SURVEY
crossdate9=find(tspan==2019);%SUBNATIONAL SURVEY
%% {Relative risks} based on littrature
RROn=[8.38,6.48]'; %RR was 7.28, 95% CI: 6.47, 8.28 for obesity by Abdullah A et al 2010 (table3 F && M)
RRSn=[1.26,1.42]';
RRFn=[1.45,1.45]';
%% estimate dm incidence rate for all risk factors 
% RRO=repmat(RROn(:,1),1,20);
% RRS=repmat(RRSn(:,1),1,20);
% RRF=repmat(RRFn(:,1),1,20);
[lambdaHn,RRO] = meshgrid(lambdaH(1,:),RROn);
[lambdaHn,RRS] = meshgrid(lambdaH(1,:),RRSn);
[lambdaHn,RRF] = meshgrid(lambdaH(1,:),RRFn);
RROS=RRO.*RRS;
RROF=RRO.*RRF;
RRSF=RRF.*RRS;
RROSF=RRO.*RRS.*RRF;
RR=[[1;1] RROn RRSn RRFn RROS(:,1) RROF(:,1) RRSF(:,1) RROSF(:,1)];
RRF(:,1:3)=0;
RRF(:,15:16)=1.48;
RRF(:,17:20)=1.30;
% {incidence rate}
lambdaO=lambdaH.*RRO;
lambdaS=lambdaH.*RRS;
lambdaF=lambdaH.*RRF;
lambdaOS=lambdaH.*RROS;
lambdaOF=lambdaH.*RROF;
lambdaSF=lambdaH.*RRSF;
lambdaOSF=lambdaH.*RROSF;

lambda=zeros(NcompTOT*intv,ag,gen);
lambdaf=[lambdaH; lambdaO; lambdaS; lambdaF; lambdaOS; lambdaOF; lambdaSF; lambdaOSF];
lambdafi=[lambdaH; lambdaO; lambdaS; lambdaF; lambdaOS; lambdaOF; lambdaSF; lambdaOSF];
lambda(:,:,1)=[lambdaf(1:2:end,:)
    zeros(8,20)
    lambdafi(1:2:end,:)
     zeros(8,20)] ;%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=[lambdaf(2:2:end,:)
    zeros(8,20)
    lambdafi(2:2:end,:)
    zeros(8,20)];%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,1,:)=0;

%%
% {Demographic parameters}
N0tot=24000;%500000;%CONSTANT POP GROWTH=1000 ELSE =10;
N0totF=0.50*N0tot;%*1.514261559;           % the total initial population
N0totM=0.50*N0tot;%*1.514261559;
N0=[N0totF,N0totM];
%% MORTALITY RATE ANALYSIS THAT ARE AGE DEPENDENT  based on litterature
muALL(:,1)=[0.003140909	0.000195455	0.000193182	0.000245455	0.000370455	0.000543182	0.000629545	0.000715909	0.001522727	0.002906818	0.005168182	0.011468182	0.024759091	0.046665909	0.069020455	0.108547727	0.168872727	0.258325	0.385511364	0.673569318]';
muALL(:,2)=[0.003610909	0.000279545	0.000272727	0.000638636	0.001218182	0.001227273	0.00155	0.001913636	0.002870455	0.004347727	0.007315909	0.013259091	0.022002273	0.045806818	0.062611364	0.100227273	0.158172727	0.244090909	0.364836364	0.619777273]';

RRm(:,1)=[5.95	5.95	5.95	5.95	5.95	5.95	5.61	5.61	3.41	3.41	2.73	2.73	2.08	2.08	 1.78	 1.78	1.78	1.78	1.78	1.78]';
RRm(:,2)=[3.70	3.70	3.70	3.70	3.70	3.70	3.30	3.30	1.95	1.95 	1.65	1.65	1.62	1.62	1.40	 1.40	1.40	1.40	1.40	1.40]';%.*(1+0.25)


%% rates for the developing different risk factors
delta(1:3,1:2)=0;
delta(4:8,1:2)=1/40;
delta(9:11,1)=1/20;
delta(9:11,2)=1/25;
delta(12:20,1)=1/7;
delta(12:20,2)=1/10;
phi=zeros(ag,gen);%0.0008; %from physical inactive to healthy

%quit smoking
epsilon=delta;%;
iota=delta;%delta;
q=delta;%delta;
tt=delta;%delta;
d=delta;%delta;
a=delta;%delta.*0;
rho=delta;%delta.*0;

%quit physical activity
theta=phi;%phi.*0;
pii=phi;%phi;
ii=phi;%phi;
p=phi;%phi;
ss=phi;%phi;
b=phi;%phi;
cc=phi;%phi;

% start smoking
nu=beta;%0;%
xi=beta;%0;%
sampi=beta;%0;%
w=beta;%beta;
n=beta;%beta;
hh=beta;%beta;
e=beta;%beta;

% start physical activity
etaa=shei;%0;%
omega=shei;%0;%
kappa=shei;%0;%
xx=shei;%shei;
ff=shei;%shei;
gg=shei;%shei;
oo=shei;%shei;

% quit obese (normal weight)
gamma=sigma;%sigma;
upsilon=sigma;%sigma;
r=sigma;%sigma;
u=sigma;%sigma;
k=sigma;%sigma;
iii=sigma;%sigma;
Leps=sigma;%sigma;
%%
% {Initial values for the variables}
% RUN FIRST TIME
x0=zeros(NcompTOT*intv,ag,gen);
for at=1:ag
        N0F(at,1)=Result(10,1)*Btot(1,1); 
        N0M(at,1)=Result(10,2)*Btot(1,2);  
end
vag=[ag,ag];
x0(1,1:ag,2)=N0M;
x0(1,1:ag,1)=N0F;

%% RUNS AFTERWARDS
% load('xinit.mat')
% x0=x(213,:)';

%%
x01=reshape(x0,numel(x0),1);
opts = odeset('Stats','on');
disp('ode45 stats:')
tic
[t,x]=ode15s(@risk_age_DM_3RFs_new,tspan,x01,opts,lambda,eta,muALL,RRm,alphas, beta, shei,sigma,scenario,Result,lambdaf,lambdafi,lambdaH,lambdaO,lambdaS,lambdaF,lambdaOS,lambdaOF,lambdaSF, lambdaOSF,intv);
toc
% {}
esp=1e-20;

%%
X=reshape(x,length(t),NcompTOT*intv,ag,gen); % X(time,stage,age, sex)
%% TOTAL POPULATION
NN = sum(x,2); %	{Population size total}
%age specific population
NGE(:,:)=sum(sum(X,3),2); % total population of each sex
NAge(:,:,:)=sum(X,2);
NAgeF=NAge(:,:,1); %	{Population size of each age group in Women}
NAgeM=NAge(:,:,2); %	{Population size of each age group in Men}
% NGE=[sum(NAgeF,2) sum(NAgeM,2)];
for agn=1:ag
    AGeDistF(:,agn)=NAgeF(:,agn)./NGE(:,1);
    AGeDistM(:,agn)=NAgeM(:,agn)./NGE(:,2);
end
TT=length(t);
S1F=NAgeF;
S1M=NAgeM;

%{Summed populations}
% S1sumHF=sum(S1F,3); %{Total susceptible 1 population size}
for at=1:ag
    S1sumF(t0-t0+1,at)=S1F(1,at);
    S1sumM(t0-t0+1,at)=S1M(1,at);
end
for tx=t0+1:2049
    for at=1:ag
        S1sumF(tx-t0+1,at)=trapz(S1F(2*(tx-t0)+1:2*(tx-t0+dt)+1,at));
         S1sumM(tx-t0+1,at)=trapz(S1M(2*(tx-t0)+1:2*(tx-t0+dt)+1,at));
    end
end

S1sumTotF=sum(S1sumF,2);
S1sumTotM=sum(S1sumM,2);
SumTOT=S1sumTotF+S1sumTotM; %total population
%% calculate population in each age group as per survey

crossdate9n=find(tx==2020);%SUBNATIONAL SURVEY
crossdate2n=find(tx==2030);%SUBNATIONAL SURVEY
crossdate3n=find(tx==2040);%SUBNATIONAL SURVEY
crossdate4n=find(tx==2049);%SUBNATIONAL SURVEY

propQatar200M=(S1sumM(2020-t0+1,1)+S1sumM(2020-t0+1,2))/S1sumTotM(2020-t0+1,1);
propQatar201M=(S1sumM(2020-t0+1,4)+S1sumM(2020-t0+1,3))/S1sumTotM(2020-t0+1,1);
propQatar202M=(S1sumM(2020-t0+1,6)+S1sumM(2020-t0+1,5))/S1sumTotM(2020-t0+1,1);
propQatar203M=(S1sumM(2020-t0+1,8)+S1sumM(2020-t0+1,7))/S1sumTotM(2020-t0+1,1);
propQatar204M=(S1sumM(2020-t0+1,10)+S1sumM(2020-t0+1,9))/S1sumTotM(2020-t0+1,1);
propQatar205M=(S1sumM(2020-t0+1,11)+S1sumM(2020-t0+1,12))/S1sumTotM(2020-t0+1,1);
propQatar206M=(S1sumM(2020-t0+1,13)+S1sumM(2020-t0+1,14))/S1sumTotM(2020-t0+1,1);
propQatar207M=(S1sumM(2020-t0+1,15)+S1sumM(2020-t0+1,16))/S1sumTotM(2020-t0+1,1);
propQatar208M=sum((S1sumM(2020-t0+1,17:20)))/S1sumTotM(2020-t0+1,1);

propQatar300M=(S1sumM(2030-t0+1,1)+S1sumM(2030-t0+1,2))/S1sumTotM(2030-t0+1,1);
propQatar301M=(S1sumM(2030-t0+1,4)+S1sumM(2030-t0+1,3))/S1sumTotM(2030-t0+1,1);
propQatar302M=(S1sumM(2030-t0+1,6)+S1sumM(2030-t0+1,5))/S1sumTotM(2030-t0+1,1);
propQatar303M=(S1sumM(2030-t0+1,8)+S1sumM(2030-t0+1,7))/S1sumTotM(2030-t0+1,1);
propQatar304M=(S1sumM(2030-t0+1,10)+S1sumM(2030-t0+1,9))/S1sumTotM(2030-t0+1,1);
propQatar305M=(S1sumM(2030-t0+1,11)+S1sumM(2030-t0+1,12))/S1sumTotM(2030-t0+1,1);
propQatar306M=(S1sumM(2030-t0+1,13)+S1sumM(2030-t0+1,14))/S1sumTotM(2030-t0+1,1);
propQatar307M=(S1sumM(2030-t0+1,15)+S1sumM(2030-t0+1,16))/S1sumTotM(2030-t0+1,1);
propQatar308M=sum((S1sumM(2030-t0+1,17:20)))/S1sumTotM(2030-t0+1,1);

propQatar400M=(S1sumM(2040-t0+1,1)+S1sumM(2040-t0+1,2))/S1sumTotM(2040-t0+1,1);
propQatar401M=(S1sumM(2040-t0+1,4)+S1sumM(2040-t0+1,3))/S1sumTotM(2040-t0+1,1);
propQatar402M=(S1sumM(2040-t0+1,6)+S1sumM(2040-t0+1,5))/S1sumTotM(2040-t0+1,1);
propQatar403M=(S1sumM(2040-t0+1,8)+S1sumM(2040-t0+1,7))/S1sumTotM(2040-t0+1,1);
propQatar404M=(S1sumM(2040-t0+1,10)+S1sumM(2040-t0+1,9))/S1sumTotM(2040-t0+1,1);
propQatar405M=(S1sumM(2040-t0+1,11)+S1sumM(2040-t0+1,12))/S1sumTotM(2040-t0+1,1);
propQatar406M=(S1sumM(2040-t0+1,13)+S1sumM(2040-t0+1,14))/S1sumTotM(2040-t0+1,1);
propQatar407M=(S1sumM(2040-t0+1,15)+S1sumM(2040-t0+1,16))/S1sumTotM(2040-t0+1,1);
propQatar408M=sum((S1sumM(2040-t0+1,17:20)))/S1sumTotM(2040-t0+1,1);

propQatar500M=(S1sumM(2049-t0+1,1)+S1sumM(2049-t0+1,2))/S1sumTotM(2049-t0+1,1);
propQatar501M=(S1sumM(2049-t0+1,4)+S1sumM(2049-t0+1,3))/S1sumTotM(2049-t0+1,1);
propQatar502M=(S1sumM(2049-t0+1,6)+S1sumM(2049-t0+1,5))/S1sumTotM(2049-t0+1,1);
propQatar503M=(S1sumM(2049-t0+1,8)+S1sumM(2049-t0+1,7))/S1sumTotM(2049-t0+1,1);
propQatar504M=(S1sumM(2049-t0+1,10)+S1sumM(2049-t0+1,9))/S1sumTotM(2049-t0+1,1);
propQatar505M=(S1sumM(2049-t0+1,11)+S1sumM(2049-t0+1,12))/S1sumTotM(2049-t0+1,1);
propQatar506M=(S1sumM(2049-t0+1,13)+S1sumM(2049-t0+1,14))/S1sumTotM(2049-t0+1,1);
propQatar507M=(S1sumM(2049-t0+1,15)+S1sumM(2049-t0+1,16))/S1sumTotM(2049-t0+1,1);
propQatar508M=sum((S1sumM(2049-t0+1,17:20)))/S1sumTotM(2049-t0+1,1);

%% ________________________________________________________________________
%__________________________________________________________________________
%% {Epidmeiologic measures}
RFAge=X;
DMALLAge=zeros(length(t),ag,gen);
PrevAgeGe=zeros(length(t),NcompTOT*intv,ag,gen);
PrevAgeGet=zeros(length(t),NcompTOT*intv,ag,gen);
PopAge=zeros(length(t),ag,gen);
  for ge=1:gen
     for zz=1:length(t)
         for kk=1:ag
DMALLAge1(zz,kk,ge)=squeeze(sum(X(zz,9:16,kk,ge))); %HDM+ODM+SDM+FDM+OSDM+OFDM+SFDM+OSFDM
DMALLAge2(zz,kk,ge)=squeeze(sum(X(zz,25:30,kk,ge))); %HDM+ODM+SDM+FDM+OSDM+OFDM+SFDM+OSFDM
DMALLAge(zz,kk,ge)=DMALLAge1(zz,kk,ge)+DMALLAge2(zz,kk,ge);
for z=1:NcompTOT*intv
    PrevAgeGe(zz,z,kk,ge)=  X(zz,z,kk,ge)'./NGE(zz,ge);
    PrevAgeGet(zz,z,kk,ge)= X(zz,z,kk,ge)'./NN(zz,1);
end
    PopAge1(zz,kk,ge)=squeeze(sum(X(zz,1:16,kk,ge)));
    PopAge2(zz,kk,ge)=squeeze(sum(X(zz,17:32,kk,ge)));
    PopAge(zz,kk,ge)=PopAge1(zz,kk,ge)+PopAge2(zz,kk,ge);
         end
    DMGE1(zz,ge)=squeeze(sum(DMALLAge1(zz,:,ge)));
    DMNGE1(zz,ge)=squeeze(sum(PopAge1(zz,:,ge)));

    DMGE2(zz,ge)=squeeze(sum(DMALLAge2(zz,5:16,ge)));
    DMNGE2(zz,ge)=squeeze(sum(PopAge2(zz,5:16,ge)));
   
    DMGE(zz,ge)=squeeze(sum(DMALLAge(zz,:,ge)));
    DMNGE(zz,ge)=squeeze(sum(PopAge(zz,:,ge)));
    DMGE2079(zz,ge)=squeeze(sum(DMALLAge(zz,5:16,ge)));
    DMNGE2079(zz,ge)=squeeze(sum(PopAge(zz,5:16,ge)));
    
    
        end
  end

DMAGE=squeeze(sum(DMALLAge,3));
PrevDMALLall=sum(DMGE,2)./(sum(DMNGE,2));
DMPOP=sum(DMGE,2);
DMAGE2079=sum(DMGE2079,2);
POPAGE2079=sum(DMNGE2079,2);
%Obesity
AgeO1=squeeze((RFAge(:,2,:,:)+ RFAge(:,5,:,:)+ RFAge(:,6,:,:)+RFAge(:,8,:,:)+RFAge(:,10,:,:)+ RFAge(:,13,:,:)+RFAge(:,14,:,:)+RFAge(:,16,:,:))); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN
AgeO2=squeeze((RFAge(:,2+NcompTOT,:,:)+ RFAge(:,5+NcompTOT,:,:)+ RFAge(:,6+NcompTOT,:,:)+RFAge(:,8+NcompTOT,:,:)+RFAge(:,10+NcompTOT,:,:)+ RFAge(:,13+NcompTOT,:,:)+RFAge(:,14+NcompTOT,:,:)+RFAge(:,16+NcompTOT,:,:))); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN
% when reducing obesity Scenario==3
if scenario>=4 && scenario<7
AgeO2=0;
end
AgeO=AgeO1+AgeO2;

PrevAgeGeO=AgeO./PopAge*100;
PrevGeOALL(:,:)=sum(AgeO,2)./sum(PopAge,2);
O=squeeze(AgeO(:,5:16,:));
ON=squeeze(PopAge(:,5:16,:));
OALL=squeeze(sum(O,2));
ONALL=squeeze(sum(ON,2));

if scenario==6.4 || scenario==6.5 || scenario==6.9
O1=squeeze(AgeO(:,4:13,:));
ON1=squeeze(PopAge(:,4:13,:));
OALL1=squeeze(sum(O1,2));
ONALL1=squeeze(sum(ON1,2));
PrevGeOALL15_64=sum(OALL1,2)./sum(ONALL1,2);
end
PrevGeOALL20_79=sum(OALL,2)./sum(ONALL,2);
if scenario>=4 && scenario<6
    figure; plot(t(233:end), PrevGeOALL20_79(233:end))
end
if scenario==6.4 || scenario==6.5 || scenario==6.9
    figure; plot(t(233:end), PrevGeOALL15_64(233:end))
end


 %Age specific analysis
AgeS1=squeeze((RFAge(:,3,:,:)+ RFAge(:,5,:,:)+ RFAge(:,7,:,:)+RFAge(:,8,:,:)+RFAge(:,11,:,:)+ RFAge(:,13,:,:)+RFAge(:,15,:,:)+RFAge(:,16,:,:))); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN
AgeS2=squeeze((RFAge(:,3+NcompTOT,:,:)+ RFAge(:,5+NcompTOT,:,:)+ RFAge(:,7+NcompTOT,:,:)+RFAge(:,8+NcompTOT,:,:)+RFAge(:,11+NcompTOT,:,:)+ RFAge(:,13+NcompTOT,:,:)+RFAge(:,15+NcompTOT,:,:)+RFAge(:,16+NcompTOT,:,:))); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN
% % when reducing smoking Scenario==3
% AgeS2=0;
AgeS=AgeS1+AgeS2;

PrevAgeGeS=AgeS./PopAge*100;
PrevGeSALL(:,:)=sum(AgeS,2)./sum(PopAge,2);
S=squeeze(AgeS(:,5:16,:));
SN=squeeze(PopAge(:,5:16,:));
SALL=squeeze(sum(S,2));
SNALL=squeeze(sum(SN,2));
PrevGeSALL20_79=sum(SALL,2)./sum(SNALL,2);

%Physical inactivity
AgeF1=squeeze((RFAge(:,4,:,:)+ RFAge(:,6,:,:)+ RFAge(:,7,:,:)+RFAge(:,8,:,:)+RFAge(:,12,:,:)+ RFAge(:,14,:,:)+RFAge(:,15,:,:)+RFAge(:,16,:,:))); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN
AgeF2=squeeze((RFAge(:,4+NcompTOT,:,:)+ RFAge(:,6+NcompTOT,:,:)+ RFAge(:,7+NcompTOT,:,:)+RFAge(:,8+NcompTOT,:,:)+RFAge(:,12+NcompTOT,:,:)+ RFAge(:,14+NcompTOT,:,:)+RFAge(:,15+NcompTOT,:,:)+RFAge(:,16+NcompTOT,:,:))); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN
% % when reducing physical activity Scenario==3
% AgeF2=0;
AgeF=AgeF1+AgeF2;

PrevAgeGeF=  AgeF./PopAge*100; 
PrevGeFALL(:,:)=sum(AgeF,2)./sum(PopAge,2);
% FALL=sum(AgeF(:,5:16,:));%prevalence of risk factor in 15-64 year old
F=squeeze(AgeF(:,5:16,:));
FN=squeeze(PopAge(:,5:16,:));
FALL=squeeze(sum(F,2));
FNALL=squeeze(sum(FN,2));
PrevGeFALL20_79=sum(FALL,2)./sum(FNALL,2);

%Population level analysis
PrevAgeGeOALL1=squeeze((PrevAgeGet(:,2,:,:)+ PrevAgeGet(:,5,:,:)+ PrevAgeGet(:,6,:,:)+PrevAgeGet(:,8,:,:)+PrevAgeGet(:,10,:,:)+ PrevAgeGet(:,13,:,:)+PrevAgeGet(:,14,:,:)+PrevAgeGet(:,16,:,:))); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN
PrevAgeGeOALL2=squeeze((PrevAgeGet(:,2+NcompTOT,:,:)+ PrevAgeGet(:,5+NcompTOT,:,:)+ PrevAgeGet(:,6+NcompTOT,:,:)+PrevAgeGet(:,8+NcompTOT,:,:)+PrevAgeGet(:,10+NcompTOT,:,:)+ PrevAgeGet(:,13+NcompTOT,:,:)+PrevAgeGet(:,14+NcompTOT,:,:)+PrevAgeGet(:,16+NcompTOT,:,:))); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN
PrevAgeGeOALL=PrevAgeGeOALL1+PrevAgeGeOALL2;
PrevGeOALL=squeeze(sum(PrevAgeGeOALL,2));

PrevAgeGeSALL1=squeeze((PrevAgeGet(:,3,:,:)+ PrevAgeGet(:,5,:,:)+ PrevAgeGet(:,7,:,:)+PrevAgeGet(:,8,:,:)+PrevAgeGet(:,11,:,:)+ PrevAgeGet(:,13,:,:)+PrevAgeGet(:,15,:,:)+PrevAgeGet(:,16,:,:))); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN
PrevAgeGeSALL2=squeeze((PrevAgeGet(:,3,:,:)+ PrevAgeGet(:,5+NcompTOT,:,:)+ PrevAgeGet(:,7+NcompTOT,:,:)+PrevAgeGet(:,8+NcompTOT,:,:)+PrevAgeGet(:,11+NcompTOT,:,:)+ PrevAgeGet(:,13+NcompTOT,:,:)+PrevAgeGet(:,15+NcompTOT,:,:)+PrevAgeGet(:,16+NcompTOT,:,:))); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN
PrevAgeGeSALL=PrevAgeGeSALL1+PrevAgeGeSALL2;
PrevGeSALL=squeeze(sum(PrevAgeGeSALL,2));

PrevAgeGeFALL1=squeeze((PrevAgeGet(:,4,:,:)+ PrevAgeGet(:,6,:,:)+ PrevAgeGet(:,7,:,:)+PrevAgeGet(:,8,:,:)+PrevAgeGet(:,12,:,:)+ PrevAgeGet(:,14,:,:)+PrevAgeGet(:,15,:,:)+PrevAgeGet(:,16,:,:))); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN
PrevAgeGeFALL2=squeeze((PrevAgeGet(:,4+NcompTOT,:,:)+ PrevAgeGet(:,6+NcompTOT,:,:)+ PrevAgeGet(:,7+NcompTOT,:,:)+PrevAgeGet(:,8+NcompTOT,:,:)+PrevAgeGet(:,12+NcompTOT,:,:)+ PrevAgeGet(:,14+NcompTOT,:,:)+PrevAgeGet(:,15+NcompTOT,:,:)+PrevAgeGet(:,16+NcompTOT,:,:))); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN
PrevAgeGeFALL=PrevAgeGeFALL1+PrevAgeGeFALL2;
PrevGeFALL=squeeze(sum(PrevAgeGeFALL,2));

%%%%
%Total prevalence 
PrevDM=sum(sum(DMALLAge,3),2)./NN;
% Prevalence amon 15-64 year old as in Qatar population. Can be changed in
% line 227-229
PrevDMALL15_64=sum(DMGE,2)./(sum(DMNGE,2));
PrevDMALL15_641=sum(DMGE1,2)./(sum(DMNGE1,2));
PrevDMALL15_642=sum(DMGE2,2)./(sum(DMNGE2,2));

PrevDMGE=DMGE./DMNGE;

%prev of 20-79 as IDF
PrevDMALL20_79=sum(DMGE2079,2)./sum(DMNGE2079,2);
PrevDMGE2079=DMGE2079./DMNGE2079;
PrevDMAgeGE=DMALLAge./NAge;

for zz=1:length(t)
     for ge=1:gen
         % for 15-64 year old
PrevGeS(zz,ge)=sum(AgeS(zz,4:13,ge))./sum(PopAge(zz,4:13,ge));
PrevGeO(zz,ge)=sum(AgeO(zz,5:13,ge))./sum(PopAge(zz,4:13,ge));
PrevGeF(zz,ge)=sum(AgeF(zz,4:13,ge))./sum(PopAge(zz,4:13,ge));
%for 20-79 year old
PrevGeS2079(zz,ge)=sum(AgeS(zz,5:16,ge))./sum(PopAge(zz,5:16,ge));
PrevGeO2079(zz,ge)=sum(AgeO(zz,5:16,ge))./sum(PopAge(zz,5:16,ge));
PrevGeF2079(zz,ge)=sum(AgeF(zz,5:16,ge))./sum(PopAge(zz,5:16,ge));
    for rfs=1:NcompTOT*intv
       Prevge(zz,rfs,ge)=sum(PrevAgeGet(zz,rfs,5:16,ge)); 
       
       Prevge2079(zz,rfs,ge)=sum(PrevAgeGet(zz,rfs,5:16,ge));
    end
              
     end
     PrevGeOALLn(zz,1)=sum(sum(AgeO(zz,4:13,:)))./sum(sum(PopAge(zz,4:13,:)));
     PrevGeSALLn(zz,1)=sum(sum(AgeS(zz,4:13,:)))./sum(sum(PopAge(zz,4:13,:)));
     PrevGeFALLn(zz,1)=sum(sum(AgeF(zz,4:13,:)))./sum(sum(PopAge(zz,4:13,:)));

end

%for 15-64 year ald
Previntv=squeeze(sum(Prevge,3));
Prev1=Previntv(:,1:16);
Prev2=Previntv(:,17:32);
Prev=Prev1+Prev2; 

 %for 20-79 year old
Prev2079intv=squeeze(sum(Prevge2079,3));
Prev20791=Prev2079intv(:,1:16);
Prev20792=Prev2079intv(:,17:32);
Prev2079=Prev20791+Prev20792;
%%
 for ge=1:gen
        for zz=1:length(t)

   %NO INTERVENTION GROUP
   Nagetot1(zz,ge)=squeeze(sum(sum(PopAge1(zz,5:16,ge))));
   OGE1(zz,ge)=squeeze(sum(sum(AgeO1(zz,5:16,ge))));
   SGE1(zz,ge)= squeeze(sum(sum(AgeS1(zz,5:16,ge))));
   FGE1(zz,ge)=squeeze(sum(sum(AgeF1(zz,5:16,ge))));

   %INTERVENTION
   Nagetot2(zz,ge)=squeeze(sum(sum(PopAge2(zz,5:16,ge))));
%     OGE2(zz,ge)=squeeze(sum(sum(AgeO2(zz,5:16,ge))));
%     SGE2(zz,ge)= squeeze(sum(sum(AgeS2(zz,5:16,ge))));
%     FGE2(zz,ge)=squeeze(sum(sum(AgeF2(zz,5:16,ge))));
%    
   %All
   Nagetot(zz,ge)=squeeze(sum(sum(PopAge(zz,5:16,ge))));
   OGE(zz,ge)=squeeze(sum(sum(AgeO(zz,5:16,ge))));
   SGE(zz,ge)= squeeze(sum(sum(AgeS(zz,5:16,ge))));
   FGE(zz,ge)=squeeze(sum(sum(AgeF(zz,5:16,ge))));

     end
  end
PrevOALL15_64=sum(OGE,2)./sum(Nagetot,2);
PrevSALL15_64=sum(SGE,2)./sum(Nagetot,2);
PrevFALL15_64=sum(FGE,2)./sum(Nagetot,2);
PrevGeO=OGE./Nagetot;
PrevGeS=SGE./Nagetot;
PrevGeF=FGE./Nagetot;

% NO INTERVENTION GROUP
PrevSALL1= (Prev1(:,3)+ Prev1(:,5)+ Prev1(:,7)+Prev1(:,8)+Prev1(:,11)+ Prev1(:,13)+Prev1(:,15)+Prev1(:,16)); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN

PrevOALL1=(Prev1(:,2)+ Prev1(:,5)+ Prev1(:,6)+Prev1(:,8)+Prev1(:,10)+ Prev1(:,13)+Prev1(:,14)+Prev1(:,16)); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN

PrevFALL1=(Prev1(:,4)+ Prev1(:,6)+ Prev1(:,7)+Prev1(:,8)+Prev1(:,12)+ Prev1(:,14)+Prev1(:,15)+Prev1(:,16)); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN

% INTERVENTION GROUP
PrevSALL2= (Prev2(:,3)+ Prev2(:,5)+ Prev2(:,7)+Prev2(:,8)+Prev2(:,11)+ Prev2(:,13)+Prev2(:,15)+Prev2(:,16)); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN

PrevOALL2=(Prev2(:,2)+ Prev2(:,5)+ Prev2(:,6)+Prev2(:,8)+Prev2(:,10)+ Prev2(:,13)+Prev2(:,14)+Prev2(:,16)); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN

PrevFALL2=(Prev2(:,4)+ Prev2(:,6)+ Prev2(:,7)+Prev2(:,8)+Prev2(:,12)+ Prev2(:,14)+Prev2(:,15)+Prev2(:,16)); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN

%ALL
PrevOALL=PrevOALL1+PrevOALL2;
PrevSALL=PrevSALL1+PrevSALL2;
PrevFALL=PrevFALL1+PrevFALL2;

PrevSum=sum(Prev,2);%S+PrevO+PrevH+PrevOS+PrevOF+PrevSF+PrevOSF

%% INCIDENCE
%{Overall incidence}
Incidencepart=zeros(length(t),NrfTOT+1,ag,gen);
IncidencepartPtot=zeros(length(t),NrfTOT+1,ag,gen);
Incidencepartintv=zeros(length(t),NrfTOT+1,ag,gen);
IncidencepartPtotintv=zeros(length(t),NrfTOT+1,ag,gen);
% lambdaT(:,:,1)=reshape(lambdaf(:,:,1),16,20,1); 
% lambdaT(:,:,2)=reshape(lambdaf(:,:,1),16,20,1);
intvsdate=find(t==2021);
IEH=ones(8,20);
for ttt=1:length(t)
       for rfc=1:NrfTOT+1
       for kk=1:ag
           for ge=1:gen
               if ttt > intvsdate
               % if scenario==2.1 || scenario==2.2 || scenario==2.3 || scenario==2.4  || scenario==2.5 
                   if scenario==2.1; IEH=0.61; elseif scenario==2.2; IEH=0.66; elseif scenario==2.3; IEH=0.68; elseif scenario==2.4; IEH=0.74; elseif scenario==2.5;IEH=0.85;

                   % 
                   elseif scenario==3.1 || scenario==3.6;  IEH(2,:)=0.72;    IEH(5,:)=0.72;    IEH(6,:)=0.72;    IEH(8,:)=0.72; lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];

                   elseif scenario==3.2 || scenario==3.7;  IEH(2,:)=0.79;    IEH(5,:)=0.79;    IEH(6,:)=0.79;    IEH(8,:)=0.79; lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];

                   elseif scenario==3.3 || scenario==3.8;  IEH(2,8:20)=0.72; IEH(5,8:20)=0.72; IEH(6,8:20)=0.72; IEH(8,8:20)=0.72;lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];

                   elseif scenario==3.4 || scenario==3.9;  IEH(2,8:20)=0.79; IEH(5,8:20)=0.79; IEH(6,8:20)=0.79; IEH(8,8:20)=0.79;lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];

                   elseif scenario==3.5 || scenario==3.11; IEH(1:8,11:20)=0.72; lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];

                       
                   elseif scenario==4.1 || scenario==4.2; IEH(1:8,:)=0.75;     lambdafin=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
% 
                   elseif scenario==4.3 || scenario==4.4  || scenario==4.5 || scenario==4.6; IEH(4,:)=1/1.48; IEH(7,:)=1/1.48; IEH(6,:)=1/1.48; IEH(8,:)=1/1.48; % OBESITY + ANY RF (REGARDLESS OF AGE)-50%COVERAGE LSM EFFECACY=0.79  
                            lambdafin=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];

                   elseif scenario==5.1 || scenario==5.3 || scenario==5.5 || scenario==5.6; IEH(1:8,:)=1; 
                            lambdafin=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
  
                   elseif scenario==5.2 || scenario==5.4 || scenario==5.7 || scenario==5.8; IEH(1:8,:)=0.93; % 
                            lambdafin=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
              
                   elseif scenario==6.1;  IEH(1:8,:)=0.93;  lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];

                   elseif scenario==6.2 || scenario==6.6;  IEH(1:8,:)=0.90; lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];
 
                   elseif scenario==6.3 || scenario==6.7;  IEH(1:8,:)=0.87; lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];

                   elseif scenario==6.4 || scenario==6.8;  IEH(1:8,:)=1;    lambdafin=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
                       
                   elseif scenario==6.5 || scenario==6.11; IEH(1:8,:)=0.93; lambdafin=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
                   
                   elseif scenario==7; 
                      IEH(2,8:20)=0.75*0.47; IEH(5,8:20)=0.75*0.47; IEH(6,8:20)=0.75*0.47; IEH(8,8:20)=0.75*0.47;
        IEH(1:8,4:13)=0.75*0.65; IEH(1,:)=0.75*0.70;IEH(3,:)=0.75*0.70;IEH(4,:)=0.75*0.70; IEH(7,:)=0.75*0.70;% OBESITY + AGE 35+ -50%COVERAGE LSM EFFECACY=0.65                                       
                        lambdafin=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
elseif scenario==8; 
    IEH(1:8,:)=.93; IEH(2,:)=0.70*.93; IEH(5,:)=0.70*.93; IEH(6,:)=0.70*.93; IEH(8,:)=0.70*.93; 
    lambdafin=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];

                   else
                       IEH=ones(8,1); lambdafin=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];

                   end

lambda(:,:,1)=[lambdaf(1:2:end,:)
    zeros(8,20)
    lambdafin(1:2:end,:)
     zeros(8,20)] ;%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=[lambdaf(2:2:end,:)
    zeros(8,20)
    lambdafin(2:2:end,:)
    zeros(8,20)];%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,1,:)=0;
               end
           
                
Incidencepart(ttt,rfc,kk,ge)= lambda(rfc,kk,ge).*squeeze(X(ttt,rfc,kk,ge)); %lambdaH*H +  lambdaO*O + lambdaS*S +lambdaF*F + lambdaOS*OS +lambdaOF*OF + lambdaSF*SF + lambdaOSF*OSF
Incidencepartintv(ttt,rfc,kk,ge)= lambda(rfc+NcompTOT,kk,ge).*squeeze(X(ttt,rfc+NcompTOT,kk,ge)); %lambdaH*H +  lambdaO*O + lambdaS*S +lambdaF*F + lambdaOS*OS +lambdaOF*OF + lambdaSF*SF + lambdaOSF*OSF
Stot(ttt,rfc,kk,ge)=X(ttt,rfc,kk,ge)+X(ttt,rfc+NcompTOT,kk,ge);

IncidencepartPtot(ttt,rfc,kk,ge)= lambda(rfc,kk,ge).*squeeze(X(ttt,rfc,kk,ge)); %lambdaH*H +  lambdaO*O + lambdaS*S +lambdaF*F + lambdaOS*OS +lambdaOF*OF + lambdaSF*SF + lambdaOSF*OSF   
IncidencepartPtotintv(ttt,rfc,kk,ge)= lambda(rfc+NcompTOT,kk,ge).*squeeze(X(ttt,rfc+NcompTOT,kk,ge));  %lambdaH*H +  lambdaO*O + lambdaS*S +lambdaF*F + lambdaOS*OS +lambdaOF*OF + lambdaSF*SF + lambdaOSF*OSF   
%            end
           end
       end
       end
       
end

Stott=sum(sum(sum(X,4),3),2);
Ssex=squeeze(sum(sum(Stot,3),2));
Incsex=squeeze(sum(sum(Incidencepart,3),2)+sum(sum(Incidencepartintv,3),2));
IncidenceRFGE=squeeze(sum(Incidencepart,3)+sum(Incidencepartintv,3));
IncidenceRF=sum(IncidenceRFGE,3);
IncidenceTOT=sum(IncidenceRF,2);
IncidenceTOTGE(:,:)=sum(IncidenceRFGE,2);
IncRate=IncidenceTOT./Stott;
IncRatesex=Incsex./Ssex;
INcrateage=squeeze(sum(Incidencepart,2)+sum(Incidencepartintv,2))./squeeze(sum(Stot,2));
%% for onle 15-64 year old
StotAge=squeeze(sum(sum(X,4),2));
IncidenceGE=sum(Incidencepart,4)+sum(Incidencepartintv,4);
IncidenceGEsex=sum(Incidencepart,2)+sum(Incidencepartintv,2);
IncidenceRFn=squeeze(sum(IncidenceGE,2));
for zz=1:length(t)
  IncidenceRFGE15_64(zz)=sum(IncidenceRFn(zz,5:16));
  IncidenceRFGE15_64sex(zz,:)=sum(IncidenceGEsex(zz,1,5:16,:),3);
    Stot15_64(zz)=sum(StotAge(zz,5:16));
end

%%
TOTINFALL(1:length(t))=cumtrapz(tspan(1:end),IncidenceRFGE15_64(1:end));
crossdate=find(t==2021);
if scenario==0
TOTINF0(crossdate:length(t))=cumtrapz(tspan(crossdate:end),IncidenceRFGE15_64(crossdate:end));
TOTINF0ALL(crossdate:length(t))=cumtrapz(tspan(crossdate:end),IncidenceTOT(crossdate:end));
else
   TOTINFINTV(crossdate:length(t))=cumtrapz(tspan(crossdate:end),IncidenceRFGE15_64(crossdate:end));
end
% Potentia cases
IncidenceRFPGE=squeeze(sum(IncidencepartPtot,3)+sum(IncidencepartPtotintv,3));
IncidenceRFP=sum(IncidenceRFPGE,3);
IncidenceP=sum(IncidenceRFP,2);
IncidencePGE=squeeze(sum(IncidenceRFPGE,2));

%% POPULATION ATTRIBUTABLE RISK
IncidRHealthy=lambdaH;
 for zzz=1:NrfTOT+1
     for kk=1:ag
   for ge=1:gen
     FracIncidAge(:,zzz,kk,ge) = Incidencepart(:,zzz,kk,ge)./IncidenceTOTGE(:,ge)*100;
     FracIncidAgeintv(:,zzz,kk,ge) = Incidencepartintv(:,zzz,kk,ge)./IncidenceTOTGE(:,ge)*100;

ACAge1(:,zzz,kk,ge)=IncidRHealthy(ge,kk).*X(:,zzz,kk,ge).*(RR(ge,zzz)-1);
ARAge1(:,zzz,kk,ge)=IncidRHealthy(ge,kk).*X(:,zzz,kk,ge).*(RR(ge,zzz)-1)./(IncidenceTOTGE(:,ge))*100;
if RR(ge,zzz)*(1-IE)<= 1 
ACAge2(:,zzz,kk,ge)=IncidRHealthy(ge,kk).*X(:,zzz+NcompTOT,kk,ge).*(1-(RR(ge,zzz).*(1-IE)));
ARAge2(:,zzz,kk,ge)=IncidRHealthy(ge,kk).*X(:,zzz+NcompTOT,kk,ge).*(1-(RR(ge,zzz).*(1-IE)))./(IncidencePGE(:,ge))*100;

elseif RR(ge,zzz)*(1-IE)>= 1 
ACAge2(:,zzz,kk,ge)=IncidRHealthy(ge,kk).*X(:,zzz+NcompTOT,kk,ge).*(RR(ge,zzz))-IncidRHealthy(ge,kk).*X(:,zzz+NcompTOT,kk,ge).*(RR(ge,zzz).*(1-IE));
ARAge2(:,zzz,kk,ge)=IncidRHealthy(ge,kk).*X(:,zzz+NcompTOT,kk,ge).*(RR(ge,zzz)-1)./(IncidenceTOTGE(:,ge))*100;
% if zzz >=8
%    ACAgeP(:,zzz,kk,ge)=IncidRHealthy(ge,kk).*X(:,zzz,kk,ge).*(RR(ge,zzz-3))-IncidRHealthy(ge,kk).*X(:,zzz,kk,ge).*(RR(ge,zzz));
end
    end
    FracIncid(:,zzz,ge) = IncidenceRFGE(:,zzz,ge)./IncidenceTOTGE(:,ge)*100;
        AR=sum(ARAge1,3);
          
     end
 end
AC(:,:,:)=sum(ACAge1,3);
ACGE(:,:)=sum(AC,3);
%% proportional risk to all RFs
 for kk=5:16 
     for ge=1:gen
PAROALLAge2ge(:,kk,ge)=(ACAge1(:,2,kk,ge)+ACAge1(:,5,kk,ge).*((RR(ge,2)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)))+ACAge1(:,6,kk,ge).*((RR(ge,2)-1)/((RR(ge,2)-1)+(RR(ge,4)-1)))+ACAge1(:,8,kk,ge).*((RR(ge,2)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)+(RR(ge,4)-1)))+ACAge2(:,2,kk,ge)+ACAge2(:,5,kk,ge)+ACAge2(:,6,kk,ge)+ACAge2(:,8,kk,ge))./(IncidenceTOTGE(:,ge))*100;%(IncidRHealthy(1,kk).*X(:,2,kk).*(RR(1,2)-1)+IncidRHealthy(1,kk).*X(:,5,kk).*(RR(1,5)-1).*((RR(1,2)-1)/((RR(1,2)-1)+(RR(1,3)-1)))+IncidRHealthy(1,kk).*X(:,6,kk).*(RR(1,6)-1).*((RR(1,2)-1)./((RR(1,2)-1)+(RR(1,4)-1)))+IncidRHealthy(1,kk).*X(:,8,kk).*(RR(1,8)-1).*((RR(1,2)-1)./((RR(1,2)-1)+(RR(1,3)-1)+(RR(1,4)-1))))./(Incidence(:,1))*100;
PARSALLAge2ge(:,kk,ge)=(ACAge1(:,3,kk,ge)+ACAge1(:,5,kk,ge).*((RR(ge,3)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)))+ACAge1(:,7,kk,ge).*((RR(ge,3)-1)/((RR(ge,4)-1)+(RR(ge,3)-1)))+ACAge1(:,8,kk,ge).*((RR(ge,3)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)+(RR(ge,4)-1)))+ACAge2(:,3,kk,ge)+ACAge2(:,5,kk,ge)+ACAge2(:,7,kk,ge)+ACAge2(:,8,kk,ge))./(IncidenceTOTGE(:,ge))*100;%(IncidRHealthy(1,kk).*X(:,3,kk).*(RR(1,3)-1)+IncidRHealthy(1,kk).*X(:,5,kk).*(RR(1,5)-1).*((RR(1,3)-1)/((RR(1,2)-1)+(RR(1,3)-1)))+IncidRHealthy(1,kk).*X(:,7,kk).*(RR(1,7)-1).*((RR(1,3)-1)./((RR(1,3)-1)+(RR(1,4)-1)))+IncidRHealthy(1,kk).*X(:,8,kk).*(RR(1,8)-1).*((RR(1,3)-1)./((RR(1,2)-1)+(RR(1,3)-1)+(RR(1,4)-1))))./(Incidence(:,1))*100;
PARFALLAge2ge(:,kk,ge)=((ACAge1(:,4,kk,ge)+ACAge1(:,6,kk,ge).*((RR(ge,4)-1)/((RR(ge,2)-1)+(RR(ge,4)-1)))+ACAge1(:,7,kk,ge).*((RR(ge,4)-1)/((RR(ge,4)-1)+(RR(ge,3)-1)))+ACAge1(:,8,kk,ge).*((RR(ge,4)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)+(RR(ge,4)-1))))+ACAge2(:,4,kk,ge)+ACAge2(:,7,kk,ge)+ACAge2(:,6,kk,ge)+ACAge2(:,8,kk,ge))./(IncidenceTOTGE(:,ge))*100;%(IncidRHealthy(1,kk).*X(:,4,kk).*(RR(1,4)-1)+IncidRHealthy(1,kk).*X(:,6,kk).*(RR(1,6)-1).*((RR(1,4)-1)/((RR(1,2)-1)+(RR(1,4)-1)))+IncidRHealthy(1,kk).*X(:,7,kk).*(RR(1,7)-1).*((RR(1,4)-1)./((RR(1,3)-1)+(RR(1,4)-1)))+IncidRHealthy(1,kk).*X(:,8,kk).*(RR(1,8)-1).*((RR(1,4)-1)./((RR(1,2)-1)+(RR(1,3)-1)+(RR(1,4)-1))))./(Incidence(:,1))*100;
 PPF(:,kk,ge)=(ACAge2(:,1,kk,ge)+ACAge2(:,2,kk,ge)+ACAge2(:,3,kk,ge)+ACAge2(:,4,kk,ge)+ACAge2(:,5,kk,ge)+ACAge2(:,7,kk,ge)+ACAge2(:,6,kk,ge)+ACAge2(:,8,kk,ge))./(IncidencePGE(:,ge))*100;

 PAROALLAge2(:,kk,ge)=(ACAge1(:,2,kk,ge)+ACAge1(:,5,kk,ge).*((RR(ge,2)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)))+ACAge1(:,6,kk,ge).*((RR(ge,2)-1)/((RR(ge,2))))+ACAge1(:,8,kk,ge).*((RR(ge,2)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)))+ACAge2(:,2,kk,ge)+ACAge2(:,5,kk,ge)+ACAge2(:,6,kk,ge)+ACAge2(:,8,kk,ge))./(IncidenceTOT(:,1))*100;%(IncidRHealthy(1,kk).*X(:,2,kk).*(RR(1,2)-1)+IncidRHealthy(1,kk).*X(:,5,kk).*(RR(1,5)-1).*((RR(1,2)-1)/((RR(1,2)-1)+(RR(1,3)-1)))+IncidRHealthy(1,kk).*X(:,6,kk).*(RR(1,6)-1).*((RR(1,2)-1)./((RR(1,2)-1)+(RR(1,4)-1)))+IncidRHealthy(1,kk).*X(:,8,kk).*(RR(1,8)-1).*((RR(1,2)-1)./((RR(1,2)-1)+(RR(1,3)-1)+(RR(1,4)-1))))./(Incidence(:,1))*100;
PARSALLAge2(:,kk,ge)=(ACAge1(:,3,kk,ge)+ACAge1(:,5,kk,ge).*((RR(ge,3)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)))+ACAge1(:,7,kk,ge).*((RR(ge,3)-1)/((RR(ge,3))))+ACAge1(:,8,kk,ge).*((RR(ge,3)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)))+ACAge2(:,3,kk,ge)+ACAge2(:,5,kk,ge)+ACAge2(:,7,kk,ge)+ACAge2(:,8,kk,ge))./(IncidenceTOT(:,1))*100;%(IncidRHealthy(1,kk).*X(:,3,kk).*(RR(1,3)-1)+IncidRHealthy(1,kk).*X(:,5,kk).*(RR(1,5)-1).*((RR(1,3)-1)/((RR(1,2)-1)+(RR(1,3)-1)))+IncidRHealthy(1,kk).*X(:,7,kk).*(RR(1,7)-1).*((RR(1,3)-1)./((RR(1,3)-1)+(RR(1,4)-1)))+IncidRHealthy(1,kk).*X(:,8,kk).*(RR(1,8)-1).*((RR(1,3)-1)./((RR(1,2)-1)+(RR(1,3)-1)+(RR(1,4)-1))))./(Incidence(:,1))*100;
PARFALLAge2(:,kk,ge)=((ACAge1(:,4,kk,ge)+ACAge1(:,6,kk,ge).*((RR(ge,4)-1)/((RR(ge,2)-1)+(RR(ge,4)-1)))+ACAge1(:,7,kk,ge).*((RR(ge,4)-1)/((RR(ge,3)-1)+(RR(ge,4)-1)))+ACAge1(:,8,kk,ge).*((RR(ge,4)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)+(RR(ge,4)-1))))+ACAge2(:,4,kk,ge)+ACAge2(:,7,kk,ge)+ACAge2(:,6,kk,ge)+ACAge2(:,8,kk,ge))./(IncidenceTOT(:,1))*100;%(IncidRHealthy(1,kk).*X(:,4,kk).*(RR(1,4)-1)+IncidRHealthy(1,kk).*X(:,6,kk).*(RR(1,6)-1).*((RR(1,4)-1)/((RR(1,2)-1)+(RR(1,4)-1)))+IncidRHealthy(1,kk).*X(:,7,kk).*(RR(1,7)-1).*((RR(1,4)-1)./((RR(1,3)-1)+(RR(1,4)-1)))+IncidRHealthy(1,kk).*X(:,8,kk).*(RR(1,8)-1).*((RR(1,4)-1)./((RR(1,2)-1)+(RR(1,3)-1)+(RR(1,4)-1))))./(Incidence(:,1))*100;
 PPFAge(:,kk,ge)=(ACAge2(:,1,kk,ge)+ACAge2(:,2,kk,ge)+ACAge2(:,3,kk,ge)+ACAge2(:,4,kk,ge)+ACAge2(:,5,kk,ge)+ACAge2(:,7,kk,ge)+ACAge2(:,6,kk,ge)+ACAge2(:,8,kk,ge))./(IncidenceP(:,1))*100;

% PAROALLAge2ex(:,kk,ge)=(ACAge(:,2,kk,ge)+ACAge(:,5,kk,ge).*((RR(ge,2)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)))+ACAge(:,6,kk,ge).*((RR(ge,2)-1)/((RR(ge,2))))+ACAge(:,8,kk,ge).*((RR(ge,2)-1)/((RR(ge,2)-1)+(RR(ge,3)-1))));
% PARSALLAge2ex(:,kk,ge)=(ACAge(:,3,kk,ge)+ACAge(:,5,kk,ge).*((RR(ge,3)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)))+ACAge(:,7,kk,ge).*((RR(ge,3)-1)/((RR(ge,3))))+ACAge(:,8,kk,ge).*((RR(ge,3)-1)/((RR(ge,2)-1)+(RR(ge,3)-1))));
% PARFALLAge2ex(:,kk,ge)=((ACAge(:,4,kk,ge)+ACAge(:,6,kk,ge).*((RR(ge,4)-1)/((RR(ge,2)-1)+(RR(ge,4)-1)))+ACAge(:,7,kk,ge).*((RR(ge,4)-1)/((RR(ge,3)-1)+(RR(ge,4)-1)))+ACAge(:,8,kk,ge).*((RR(ge,4)-1)/((RR(ge,2)-1)+(RR(ge,3)-1)+(RR(ge,4)-1)))));    

     end
 end
%denominator is incidence in the specific sex (more correct)
PAROALL2ge(:,:)=sum(PAROALLAge2ge,2);
PARSALL2ge(:,:)=sum(PARSALLAge2ge,2);
PARFALL2ge(:,:)=sum(PARFALLAge2ge,2);

PAROALL2geN=squeeze(PAROALLAge2ge(:, 5:16,:));
PARSALL2geN=squeeze(PARSALLAge2ge(:, 5:16,:));
PARFALL2geN=squeeze(PARFALLAge2ge(:, 5:16,:));


PAROALL2ge20_79(:,:)=squeeze(sum(PAROALL2geN,2));
PARSALL2ge20_79(:,:)=squeeze(sum(PARSALL2geN,2));
PARFALL2ge20_79(:,:)=squeeze(sum(PARFALL2geN,2));

%denominator is incidence in the total pop (sum age and sex)
PAROALL2=sum(sum(PAROALLAge2,3),2);
PARSALL2=sum(sum(PARSALLAge2,3),2);
PARFALL2=sum(sum(PARFALLAge2,3),2);
PARALLALL2=PAROALL2+PARSALL2+PARFALL2;

PAROALL2N=(PAROALLAge2(:,5:16,:));
PARSALL2N=(PARSALLAge2(:,5:16,:));
PARFALL2N=(PARFALLAge2(:,5:16,:));


PAROALL220_79=sum(sum(PAROALL2N,3),2);
PARSALL220_79=sum(sum(PARSALL2N,3),2);
PARFALL220_79=sum(sum(PARFALL2N,3),2);
% PARALLALL2=PAROALL2+PARSALL2+PARFALL2;
%% =====================================================================================

%%
PREVF=(DMALLAge(crossdate,:,1))./(PopAge(crossdate,:,1));
PREVM=DMALLAge(crossdate,:,2)./PopAge(crossdate,:,2);

for df=1:length(t); origind(1,df)=sum(sum(sum(X(df,1:16,4:13,:)))); cov(1,df)=sum(sum(sum(X(df,17:32,4:13,:)))); end
coverage=cov./(origind+cov);
figure
plot(t,coverage.*100)

runs=length(t);
DMmort=zeros(runs,20,2,32);
for zx=1:runs
     xn(:,:)=x(zx,:);
    t1=tspan(zx);
%  
 [dxL,alpha,dmmort]=risk_age_DM_3RFs_new(t1,xn,lambda,eta,muALL,RRm,alphas, beta, shei,sigma,scenario,Result,lambdaf,lambdafi,lambdaH,lambdaO,lambdaS,lambdaF,lambdaOS,lambdaOF,lambdaSF, lambdaOSF,intv); 
%     alphasn(zx,:,:)=alpha;
    DMmort(zx,:,:,:)=dmmort;
%     Mort(zx,:,:)=mort;
end%    
% end
DMmortgr=sum(DMmort,4);
DMmortAll=sum(DMmortgr,3);
DmMortALL15_64=sum(DMmortAll,2);
TOTmort0(235:length(t))=cumtrapz(tspan(235:end),DmMortALL15_64(235:end,1));

out=[PrevDMALL15_64,PrevDMGE,IncidenceRFGE15_64',IncidenceRFGE15_64sex,IncRate,TOTINFINTV',TOTmort0'];
SHOWRES(1,:)=out(243,:);SHOWRES(2,:)=out(253,:);SHOWRES(3,:)=out(273,:);SHOWRES(4,:)=out(293,:);

