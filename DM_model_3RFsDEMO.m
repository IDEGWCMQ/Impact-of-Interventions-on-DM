function [propQatar101F,propQatar102F,propQatar103F,propQatar104F,propQatar105F,propQatar106F,...  
    propQatar121F,propQatar125F,...
    propQatar131F,propQatar132F,propQatar133F,propQatar134F,propQatar135F,...
    propQatar151F,propQatar152F,propQatar155F,...
    propQatar190F,propQatar191F,propQatar192F,propQatar193F,propQatar194F,propQatar195F,propQatar196F,propQatar197F,propQatar198F,...
    propQatar101M,propQatar102M,propQatar103M,propQatar104M,propQatar105M,propQatar106M,...  
    propQatar121M,propQatar125M,...
    propQatar131M,propQatar132M,propQatar133M,propQatar134M,propQatar135M,...
    propQatar151M,propQatar152M,propQatar155M,...
    propQatar190M,propQatar191M,propQatar192M,propQatar193M,propQatar194M,propQatar195M,propQatar196M,propQatar197M,propQatar198M,...
    SumTOT,S1sumTotF,S1sumTotM]=DM_model_3RFsDEMO(Result,fittingparam1,Btot)
%DM MODEL  September 07, 2015
%SUSANNE AWAD
global ag gen NcompTOT N0 rf NrfTOT delta phi epsilon iota q tt d a rho theta pii ii p ss b cc nu xi sampi w n hh e etaa omega kappa xx ff gg oo gamma upsilon r u k iii Leps chi psi Omg m jj v l 
%% SEX 1: FEMALE, 2: MALE
gen=2;
%% Fitting paramETER
param1s=fittingparam1(1:20*2*3,1);
param1=reshape(param1s,20,2,3);
% SMOKING RATE
beta=param1(:,:,1);
% beta(1:3,:)=0;
% PHIA RATE
shei=param1(:,:,2);
sigma=param1(:,:,3);
% DM INCIDENCE RATE
param3s=fittingparam1(20*2*3+1:(20*2*3)+(6),1,1);
param3=reshape(param3s,3,2,1);
lambdaHh=param3;
% ASSUMING AGE-SPECIFIC RATES FOLLOW A GAUSSIAN DISTRIBUTION 
lambdaH=(DMincFuncGaussian(ag,gen,lambdaHh(1,:),lambdaHh(2,:),lambdaHh(3,:)));
% lambdaH(:,16:20)=0;
% OBESITY RATE
param2s=fittingparam1((20*2*3)+(6)+1:(20*2*3+6)+(3*2),1);
param2=reshape(param2s,3,2,1);
alphas=param2;

%% DEMONGRAFIC PARAMETERS
%% IN CASE OF INTERVENTIONS SCENARIO NO. 1: NO INTERVENTION
scenario=1;
%% ========================================================================
%%INPUT PARAMETER
%%=========================================================================
%% Time scale
t0=1904;         %Start time
tf=2050;         %Stop time
dt=0.5;         %time interval
tspan=t0:dt:tf;  %timespan to use in the ode45
%% CROSS SECTIONAL YEAR FOR INVESTIGATION (BASED ON SURVEY YEAR)
crossdate1=find(tspan==2017);
crossdate2=find(tspan==2009);
crossdate9=find(tspan==2040);

%% AGR STRATIFICATION
age_gr=5;       % 5-year age stratification of the age groups
ag=100/age_gr;   % Number of age groups
Dsys=100;       % The maximum age included
eta(1:ag,1)=ag/Dsys;       % ageing rate
eta(ag,1)=0;
%% {Relative risks OF RISK FACTORS}
RROn=[8.38,6.48]'; %RR was 7.28, 95% CI: 6.47, 8.28 for obesity by Abdullah A et al 2010 (table3 F && M)
RRSn=[1.33,1.42]';
RRFn=[1.45,1.45]';
%% INCIDENCE RATE
[lambdaHn,RRO] = meshgrid(lambdaH(1,:),RROn);   %FOR INDIVIDUALS THAT ARE OBESE
[lambdaHn,RRS] = meshgrid(lambdaH(1,:),RRSn);   %FOR INDIVIDUALS THAT ARE SMOKERS
[lambdaHn,RRF] = meshgrid(lambdaH(1,:),RRFn);   %FOR INDIVIDUALS THAT ARE PHIA
RROS=RRO.*RRS;              %FOR INDIVIDUALS THAT ARE OBESE AMD SMOKERS
RROF=RRO.*RRF;              %FOR INDIVIDUALS THAT ARE OBESE AND PHIA
RRSF=RRF.*RRS;              %FOR INDIVIDUALS THAT ARE SMOKERS AND PHIA
RROSF=RRO.*RRS.*RRF;        %FOR INDIVIDUALS THAT ARE OBESE SMOKERS AND PHIA
RR=[[1;1] RROn RRSn RRFn RROS(:,1) RROF(:,1) RRSF(:,1) RROSF(:,1)]; % ALL
RRF(:,1:3)=0;
RRF(:,15:16)=1.48;
RRF(:,17:20)=1.30;
% RISK FACTOR SPECIFIC DM INCIDENCE RATE
lambdaO=lambdaH.*RRO;
lambdaS=lambdaH.*RRS;
lambdaF=lambdaH.*RRF;
lambdaOS=lambdaH.*RROS;
lambdaOF=lambdaH.*RROF;
lambdaSF=lambdaH.*RRSF;
lambdaOSF=lambdaH.*RROSF;

lambdaf=[lambdaH; lambdaO; lambdaS; lambdaF; lambdaOS; lambdaOF; lambdaSF; lambdaOSF];
lambda(:,:,1)=lambdaf(1:2:end,:);%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=lambdaf(2:2:end,:);%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);

%% {Demographic parameters}
N0tot=24000;%500000;%CONSTANT POP GROWTH=1000 ELSE =10;
N0totF=0.50*N0tot;%*1.514261559;           % the total initial population
N0totM=0.50*N0tot;%*1.514261559;
N0=[N0totF,N0totM];
% MORTALITY RATE ANALYSIS THAT ARE AGE DEPENDENT  (UNCOMENT BELOW IF NEEDED)
muALL(:,1)=[0.003140909	0.000195455	0.000193182	0.000245455	0.000370455	0.000543182	0.000629545	0.000715909	0.001522727	0.002906818	0.005168182	0.011468182	0.024759091	0.046665909	0.069020455	0.108547727	0.168872727	0.258325	0.385511364	0.673569318]';
muALL(:,2)=[0.003610909	0.000279545	0.000272727	0.000638636	0.001218182	0.001227273	0.00155	0.001913636	0.002870455	0.004347727	0.007315909	0.013259091	0.022002273	0.045806818	0.062611364	0.100227273	0.158172727	0.244090909	0.364836364	0.619777273]';

%I had to cahne the last RR's for the prevalence of DM and RF to decline
RRDECODA(:,1)=[5.95	5.95	5.95	5.95	5.95	5.95	5.61	5.61	3.41	3.41	2.73	2.73	2.08	2.08	 1.78	 1.78	1.78	1.78	1.78	1.78]';
RRDECODA(:,2)=[3.70	3.70	3.70	3.70	3.70	3.70	3.30	3.30	1.95	1.95 	1.65	1.65	1.62	1.62	1.40	 1.40	1.40	1.40	1.40	1.40]';

%POPULATION GROWTH PARAMETERS (birth rate among qataris Only)
%% OTHERS
rf=3;         % Number of risk factors (from healthy to obese, smoker and physical inactivity)
NrfTOT=2^(rf)-1; % total number of risk factor with overlapping
NcompTOT=NrfTOT*2+2; %total number of compartment that the population is stratified by (including healthy)

%% RATES OF DEVELOPING OR REVERSE RISK FACTORS 
% start smoking
nu=beta;%0;%
xi=beta;%0;%
sampi=beta;%0;%
w=beta;%beta;
n=beta;%beta;
hh=beta;%beta;
e=beta;%beta;

% from SMOKER to health
delta(1:3,1:2)=0;
delta(4:8,1:2)=1/40;
delta(9:11,1)=1/20;
delta(9:11,2)=1/25;
delta(12:20,1)=1/7;
delta(12:20,2)=1/10;
%quit smoking
epsilon=delta;%;
iota=delta;%delta;
q=delta;%delta;
tt=delta;%delta;
d=delta;%delta;
a=delta;%delta.*0;
rho=delta;%delta.*0;

% start physical activity
etaa=shei;%0;%
omega=shei;%0;%
kappa=shei;%0;%
xx=shei;%shei;
ff=shei;%shei;
gg=shei;%shei;
oo=shei;%shei;

%from physical inactive to healthy
phi=zeros(ag,gen); 
%quit physical activity
theta=phi;%phi.*0;
pii=phi;%phi;
ii=phi;%phi;
p=phi;%phi;
ss=phi;%phi;
b=phi;%phi;
cc=phi;%phi;

% quit obese (normal weight)
gamma=sigma;%sigma;
upsilon=sigma;%sigma;
r=sigma;%sigma;
u=sigma;%sigma;
k=sigma;%sigma;
iii=sigma;%sigma;
Leps=sigma;%sigma;
%% ========================================================================
%%ODE CALCULATION
%%=========================================================================
%% {Initial values for the variables}
x0=zeros(NcompTOT,ag,gen);
for at=1:ag
        N0F(at,1)=Result(10,1)*Btot(1,1); 
        N0M(at,1)=Result(10,2)*Btot(1,2); 
%         N0(a,2)=0.6883*B(1,a+1);
%         N0(a,3)=0.2304*B(1,a+1);
end
vag=[ag,ag];
x0(1,1:ag,2)=N0M;
x0(1,1:ag,1)=N0F;
x01=reshape(x0,numel(x0),1);
tic
[t,x]=ode15s(@risk_age_DM_3RFs_new,tspan,x01,[],lambda,eta,muALL,RRDECODA,alphas, beta, shei,sigma, scenario, Result);
toc
esp=1e-20;
%%
X=reshape(x,length(t),NcompTOT,ag,gen);

%% TOTAL POPULATION
%age specific population
NGE(:,:)=sum(sum(X,3),2);
NAge(:,:,:)=squeeze(sum(X,2));
NAgeF=squeeze(NAge(:,:,1));
NAgeM=squeeze(NAge(:,:,2));
% NGE=[sum(NAgeF,2) sum(NAgeM,2)];
for agn=1:ag
    AGeDistF(:,agn)=NAgeF(:,agn)./NGE(:,1);
    AGeDistM(:,agn)=NAgeM(:,agn)./NGE(:,2);
end
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
SumTOT=S1sumTotF+S1sumTotM;

propQatar101F=(S1sumF(2010-t0+1,4)+S1sumF(2010-t0+1,5))/S1sumTotF(2010-t0+1,1);
propQatar102F=(S1sumF(2010-t0+1,6)+S1sumF(2010-t0+1,7))/S1sumTotF(2010-t0+1,1);
propQatar103F=(S1sumF(2010-t0+1,8)+S1sumF(2010-t0+1,9))/S1sumTotF(2010-t0+1,1);
propQatar104F=(S1sumF(2010-t0+1,10)+S1sumF(2010-t0+1,11))/S1sumTotF(2010-t0+1,1);
propQatar105F=(S1sumF(2010-t0+1,12)+S1sumF(2010-t0+1,13))/S1sumTotF(2010-t0+1,1);
propQatar106F=sum((S1sumF(2010-t0+1,14:20)))/S1sumTotF(2010-t0+1,1);
%2012
for age=1:10
propQatar121F(age)=(S1sumF(2012-t0+1,age+3))/S1sumTotF(2012-t0+1,1);
end
propQatar125F=sum((S1sumF(2012-t0+1,14:20)))/S1sumTotF(2012-t0+1,1);

%2013
propQatar131F=(S1sumF(2013-t0+1,4)+S1sumF(2013-t0+1,5))/S1sumTotF(2013-t0+1,1);
propQatar132F=(S1sumF(2013-t0+1,6)+S1sumF(2013-t0+1,7))/S1sumTotF(2013-t0+1,1);
propQatar133F=(S1sumF(2013-t0+1,8)+S1sumF(2013-t0+1,9))/S1sumTotF(2013-t0+1,1);
propQatar134F=(S1sumF(2013-t0+1,10)+S1sumF(2013-t0+1,11))/S1sumTotF(2013-t0+1,1);
propQatar135F=sum((S1sumF(2013-t0+1,12:20)))/S1sumTotF(2013-t0+1,1);
%2015
propQatar151F=sum((S1sumF(2015-t0+1,1:3)))/S1sumTotF(2015-t0+1,1);

for age=1:12
propQatar152F(age)=(S1sumF(2015-t0+1,age+3))/S1sumTotF(2015-t0+1,1);
end
propQatar155F=sum((S1sumF(2015-t0+1,16:20)))/S1sumTotF(2015-t0+1,1);
%2019
propQatar190F=(S1sumF(2019-t0+1,1)+S1sumF(2019-t0+1,2))/S1sumTotF(2019-t0+1,1);
propQatar191F=(S1sumF(2019-t0+1,4)+S1sumF(2019-t0+1,3))/S1sumTotF(2019-t0+1,1);
propQatar192F=(S1sumF(2019-t0+1,6)+S1sumF(2019-t0+1,5))/S1sumTotF(2019-t0+1,1);
propQatar193F=(S1sumF(2019-t0+1,8)+S1sumF(2019-t0+1,7))/S1sumTotF(2019-t0+1,1);
propQatar194F=(S1sumF(2019-t0+1,10)+S1sumF(2019-t0+1,9))/S1sumTotF(2019-t0+1,1);
propQatar195F=(S1sumF(2019-t0+1,11)+S1sumF(2019-t0+1,12))/S1sumTotF(2019-t0+1,1);
propQatar196F=(S1sumF(2019-t0+1,13)+S1sumF(2019-t0+1,14))/S1sumTotF(2019-t0+1,1);
propQatar197F=(S1sumF(2019-t0+1,15)+S1sumF(2019-t0+1,16))/S1sumTotF(2019-t0+1,1);
propQatar198F=sum((S1sumF(2019-t0+1,17:20)))/S1sumTotF(2019-t0+1,1);


%%
propQatar101M=(S1sumM(2010-t0+1,4)+S1sumM(2010-t0+1,5))/S1sumTotM(2010-t0+1,1);
propQatar102M=(S1sumM(2010-t0+1,6)+S1sumM(2010-t0+1,7))/S1sumTotM(2010-t0+1,1);
propQatar103M=(S1sumM(2010-t0+1,8)+S1sumM(2010-t0+1,9))/S1sumTotM(2010-t0+1,1);
propQatar104M=(S1sumM(2010-t0+1,10)+S1sumM(2010-t0+1,11))/S1sumTotM(2010-t0+1,1);
propQatar105M=(S1sumM(2010-t0+1,12)+S1sumM(2010-t0+1,13))/S1sumTotM(2010-t0+1,1);
propQatar106M=sum((S1sumM(2010-t0+1,14:20)))/S1sumTotM(2010-t0+1,1);
%2012
for age=1:10
propQatar121M(age)=(S1sumM(2012-t0+1,age+3))/S1sumTotM(2012-t0+1,1);
end
propQatar125M=sum((S1sumM(2012-t0+1,14:20)))/S1sumTotM(2012-t0+1,1);

%2013
propQatar131M=(S1sumM(2013-t0+1,4)+S1sumM(2013-t0+1,5))/S1sumTotM(2013-t0+1,1);
propQatar132M=(S1sumM(2013-t0+1,6)+S1sumM(2013-t0+1,7))/S1sumTotM(2013-t0+1,1);
propQatar133M=(S1sumM(2013-t0+1,8)+S1sumM(2013-t0+1,9))/S1sumTotM(2013-t0+1,1);
propQatar134M=(S1sumM(2013-t0+1,10)+S1sumM(2013-t0+1,11))/S1sumTotM(2013-t0+1,1);
propQatar135M=sum((S1sumM(2013-t0+1,12:20)))/S1sumTotM(2013-t0+1,1);
%2015
propQatar151M=sum((S1sumM(2015-t0+1,1:3)))/S1sumTotM(2015-t0+1,1);

for age=1:12
propQatar152M(age)=(S1sumM(2015-t0+1,age+3))/S1sumTotM(2015-t0+1,1);
end
propQatar155M=sum((S1sumM(2015-t0+1,16:20)))/S1sumTotM(2015-t0+1,1);
%2019
propQatar190M=(S1sumM(2019-t0+1,1)+S1sumM(2019-t0+1,2))/S1sumTotM(2019-t0+1,1);
propQatar191M=(S1sumM(2019-t0+1,4)+S1sumM(2019-t0+1,3))/S1sumTotM(2019-t0+1,1);
propQatar192M=(S1sumM(2019-t0+1,6)+S1sumM(2019-t0+1,5))/S1sumTotM(2019-t0+1,1);
propQatar193M=(S1sumM(2019-t0+1,8)+S1sumM(2019-t0+1,7))/S1sumTotM(2019-t0+1,1);
propQatar194M=(S1sumM(2019-t0+1,10)+S1sumM(2019-t0+1,9))/S1sumTotM(2019-t0+1,1);
propQatar195M=(S1sumM(2019-t0+1,11)+S1sumM(2019-t0+1,12))/S1sumTotM(2019-t0+1,1);
propQatar196M=(S1sumM(2019-t0+1,13)+S1sumM(2019-t0+1,14))/S1sumTotM(2019-t0+1,1);
propQatar197M=(S1sumM(2019-t0+1,15)+S1sumM(2019-t0+1,16))/S1sumTotM(2019-t0+1,1);
propQatar198M=sum((S1sumM(2019-t0+1,17:20)))/S1sumTotM(2019-t0+1,1);

% PryM=NAgeM(crossdate9,:)./(sum(NAgeM(crossdate9,:))+sum(NAgeF(crossdate9,:)))*1000000;
end