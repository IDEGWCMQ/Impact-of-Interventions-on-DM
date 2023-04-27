function [outy, S1sumF, S1sumM]=DM_model_3RFsDMQatar(fittingparam1, Result, BFf, BMm,Btot)
%DM MODEL  September 07, 2015 
%% update 2021
%SUSANNE AWAD
global ag gen NcompTOT N0 rf NrfTOT delta phi epsilon iota q tt d a rho theta pii ii p ss b cc nu xi sampi w n hh e etaa omega kappa xx ff gg oo gamma upsilon r u k iii Leps 

%%
age_gr=5;       % 5-year age stratification of the age groups
ag=100/age_gr;   % Number of age groups
Dsys=100;       % The maximum age included
eta(1:ag,1)=ag/Dsys;       % ageing rate
eta(ag,1)=0;
gen=2;
%% %%fittingparam
param1s=fittingparam1(1:20*2*3,1);
param1=reshape(param1s,20,2,3);
beta=param1(:,:,1);
shei=param1(:,:,2);
sigma=param1(:,:,3);

param3s=fittingparam1(20*2*3+1:(20*2*3)+(6),1,1);
param3=reshape(param3s,3,2,1);
lambdaHh=param3;

lambdaH=DMincFuncGaussian(ag,gen,lambdaHh(1,:),lambdaHh(2,:),lambdaHh(3,:));
 lambdaH(:,1:3)=0;

param2s=fittingparam1((20*2*3)+(6)+1:(20*2*3+6)+(3*2),1);
param2=reshape(param2s,3,2,1);
alphas=param2;

scenario=1;
%% Time scale
t0=1904;         %Start time
tf=2050;         %Stop time
dt=0.5;         %time interval
tspan=t0:dt:tf;  %timespan to use in the ode45
%% SURVEY DATES
crossdate1=find(tspan==2012); %NATIONAL SURVEY
%% {Relative risks}
RROn=[8.38,6.48]'; %RR was 7.28, 95% CI: 6.47, 8.28 for obesity by Abdullah A et al 2010 (table3 F && M)
RRSn=[1.33,1.42]';
RRFn=[1.45,1.45]';
%%
RRO=repmat(RROn(:,1),1,20);
RRS=repmat(RRSn(:,1),1,20);
RRF=repmat(RRFn(:,1),1,20);
RROS=RRO.*RRS;
RROF=RRO.*RRF;
RRSF=RRF.*RRS;
RROSF=RRO.*RRS.*RRF;
RR=[[1;1] RROn RRSn RRFn RROS(:,1) RROF(:,1) RRSF(:,1) RROSF(:,1)];
RRF(:,1:3)=0;
RRF(:,15:16)=1.48;
RRF(:,17:20)=1.30;

% {Force of infection}
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

%%
% {Demographic parameters}
N0tot=24000;%500000;%CONSTANT POP GROWTH=1000 ELSE =10;
N0totF=0.50*N0tot;%*1.514261559;           % the total initial population
N0totM=0.50*N0tot;%*1.514261559;
N0=[N0totF,N0totM];
%%
% MORTALITY RATE ANALYSIS THAT ARE AGE DEPENDENT  (UNCOMENT BELOW IF
% NEEDED)
muALL(:,1)=[0.003140909	0.000195455	0.000193182	0.000245455	0.000370455	0.000543182	0.000629545	0.000715909	0.001522727	0.002906818	0.005168182	0.011468182	0.024759091	0.046665909	0.069020455	0.108547727	0.168872727	0.258325	0.385511364	0.673569318]';
muALL(:,2)=[0.003610909	0.000279545	0.000272727	0.000638636	0.001218182	0.001227273	0.00155	0.001913636	0.002870455	0.004347727	0.007315909	0.013259091	0.022002273	0.045806818	0.062611364	0.100227273	0.158172727	0.244090909	0.364836364	0.619777273]';

%I had to cahne the last RR's for the prevalence of DM and RF to decline
RRDECODA(:,1)=[5.95	5.95	5.95	5.95	5.95	5.95	5.61	5.61	3.41	3.41	2.73	2.73	2.08	2.08	 1.78	 1.78	1.78	1.78	1.78	1.78]';
RRDECODA(:,2)=[3.70	3.70	3.70	3.70	3.70	3.70	3.30	3.30	1.95	1.95 	1.65	1.65	1.62	1.62	1.40	 1.40	1.40	1.40	1.40	1.40]';
%POPULATION GROWTH PARAMETERS (birth rate among qataris Only)
% b1=0.035744039502806205472082487504; %0.02725;%UNCOMMENT NUMBER DEPENDING ON THE ANALYSIS
% %USE ONLY WITH POP GROWTH
% b6=0.095;
% b2=0.02075;
% b3=0.115025;
% b4=0.023151205;
% b5=0.03012050;

%%
% gnd=2;          % Gender; 1:Femal, 2: uncircumcised men, 3: Circumcised men
rf=3;         % Number of risk factors (from healthy to obese, smoker and physical inactivity)
NrfTOT=2^(rf)-1; % total number of risk factor with overlapping
NcompTOT=NrfTOT*2+2; %total number of compartment that the population is stratified by (including healthy)

%%
% {Rates}
% from obese to health
delta(1:3,1:2)=0;
delta(4:8,1:2)=1/40;
delta(9:11,1)=1/20;
delta(9:11,2)=1/25;
delta(12:20,1)=1/7;
delta(12:20,2)=1/10;
phi=zeros(ag,gen); %from physical inactive to healthy

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

%% {Initial values for the variables}
x0=zeros(NcompTOT,ag,gen);
for at=1:ag
      N0F(at,1)=Result(10,1)*Btot(1,1); 
        N0M(at,1)=Result(10,2)*Btot(1,2);

end
vag=[ag,ag];
x0(1,1:ag,2)=N0M;
x0(1,1:ag,1)=N0F;
x01=reshape(x0,numel(x0),1);
tic
[t,x]=ode15s(@risk_age_DM_3RFs_new,tspan,x01,[],lambda,eta,muALL,RRDECODA,alphas, beta, shei,sigma, scenario, Result);
toc
% {}
esp=1e-20;

%%
X=reshape(x,length(t),NcompTOT,ag,gen);
%% TOTAL POPULATION
NN = sum(x,2); %	{Population size of each age group}
%age specific population
NGE(:,:)=sum(sum(X,3),2);
NAge(:,:,:)=squeeze(sum(X,2));
NAgeF=squeeze(NAge(:,:,1));
NAgeM=squeeze(NAge(:,:,2));
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

% %%in 2012
%% {Epidmeiologic measures}
RFAge=X;
ll=0;
DMALLAge=zeros(length(t),ag,gen);
PrevAgeGe=zeros(length(t),NcompTOT,ag,gen);
PrevAgeGet=zeros(length(t),NcompTOT,ag,gen);
PopAge=zeros(length(t),ag,gen);
  for ge=1:gen
     for kk=1:ag
        for zz=1:length(t)
DMALLAge(zz,kk,ge)=sum(X(zz,9:16,kk,ge))'; %HDM+ODM+SDM+FDM+OSDM+OFDM+SFDM+OSFDM
for z=1:NcompTOT
    PrevAgeGe(zz,z,kk,ge)=  X(zz,z,kk,ge)'./NGE(zz,ge);
    PrevAgeGet(zz,z,kk,ge)=  X(zz,z,kk,ge)'./NN(zz,1);
end
    PopAge(zz,kk,ge)=sum(X(zz,1:16,kk,ge));
    DMGE(zz,ge)=sum(DMALLAge(zz,4:13,ge));
    DMNGE(zz,ge)=sum(NAge(zz,4:13,ge));
    DMGE2079(zz,ge)=sum(DMALLAge(zz,5:16,ge));  %USE THIS FOR DM PREVALENCE
    DMNGE2079(zz,ge)=sum(NAge(zz,5:16,ge));     %USE THIS FOR DM PREVALENCE
%     DMGE(zz,:,2)=sum(DMALLAge(zz, 4:13,2));
        end
        ll=ll+NcompTOT;
     end
  end

DMAGE(:,:)=sum(DMALLAge,3);

%% Age specific analysis FOR RISK FACTOR
%SMOKING
AgeS1=(RFAge(:,3,:,:)+ RFAge(:,5,:,:)+ RFAge(:,7,:,:)+RFAge(:,8,:,:)+RFAge(:,11,:,:)+ RFAge(:,13,:,:)+RFAge(:,15,:,:)+RFAge(:,16,:,:)); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN
AgeS(:,:,:)=AgeS1(:,1,:,:);
PrevAgeGeS=AgeS./PopAge*100;
PrevGeSALL(:,:)=sum(AgeS,2)./sum(PopAge,2);

%OBESITY
AgeO1=(RFAge(:,2,:,:)+ RFAge(:,5,:,:)+ RFAge(:,6,:,:)+RFAge(:,8,:,:)+RFAge(:,10,:,:)+ RFAge(:,13,:,:)+RFAge(:,14,:,:)+RFAge(:,16,:,:)); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN
AgeO(:,:,:)=AgeO1(:,1,:,:);
PrevAgeGeO=AgeO./PopAge*100;
PrevGeOALL(:,:)=sum(AgeO,2)./sum(PopAge,2);

%PHIA
AgeF1=(RFAge(:,4,:,:)+ RFAge(:,6,:,:)+ RFAge(:,7,:,:)+RFAge(:,8,:,:)+RFAge(:,12,:,:)+ RFAge(:,14,:,:)+RFAge(:,15,:,:)+RFAge(:,16,:,:)); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN
AgeF(:,:,:)=AgeF1(:,1,:,:);
PrevAgeGeF=  AgeF./PopAge*100; 
PrevGeFALL(:,:)=sum(AgeF,2)./sum(PopAge,2);

%prevalence of risk factor in 15-64 year old
%Population level analysis
%OBESITY
PrevAgeGeOALL1=(PrevAgeGet(:,2,:,:)+ PrevAgeGet(:,5,:,:)+ PrevAgeGet(:,6,:,:)+PrevAgeGet(:,8,:,:)+PrevAgeGet(:,10,:,:)+ PrevAgeGet(:,13,:,:)+PrevAgeGet(:,14,:,:)+PrevAgeGet(:,16,:,:)); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN
PrevAgeGeOALL(:,:,:)=PrevAgeGeOALL1(:,1,:,:);
PrevGeOALL(:,:)=sum(PrevAgeGeOALL,2);
%SMOKING
PrevAgeGeSALL1=(PrevAgeGet(:,3,:,:)+ PrevAgeGet(:,5,:,:)+ PrevAgeGet(:,7,:,:)+PrevAgeGet(:,8,:,:)+PrevAgeGet(:,11,:,:)+ PrevAgeGet(:,13,:,:)+PrevAgeGet(:,15,:,:)+PrevAgeGet(:,16,:,:)); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN
PrevAgeGeSALL(:,:,:)=PrevAgeGeSALL1(:,1,:,:);
PrevGeSALL(:,:)=sum(PrevAgeGeSALL,2);
%PHIA
PrevAgeGeFALL1=(PrevAgeGet(:,4,:,:)+ PrevAgeGet(:,6,:,:)+ PrevAgeGet(:,7,:,:)+PrevAgeGet(:,8,:,:)+PrevAgeGet(:,12,:,:)+ PrevAgeGet(:,14,:,:)+PrevAgeGet(:,15,:,:)+PrevAgeGet(:,16,:,:)); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN
PrevAgeGeFALL(:,:,:)=PrevAgeGeFALL1(:,1,:,:);
PrevGeFALL(:,:)=sum(PrevAgeGeFALL,2);

%% Total prevalence
% DM IN TOTAL POPULATION 0-99
PrevDM=sum(sum(DMALLAge,3),2)./NN;
% Prevalence amon 15-64 year old as in Qatar population. Can be changed in % line 227-229
PrevDMALL15_64=sum(DMGE,2)./sum(DMNGE,2);
PrevDMGE=DMGE./DMNGE;

%prev of 20-79 as IDF
PrevDMALL20_79=squeeze(sum(DMGE2079,2))./squeeze(sum(DMNGE2079,2));
PrevDMGE2079=DMGE2079./DMNGE2079;

PrevDMAgeGE=DMALLAge./NAge;
for zz=1:length(t)
     for ge=1:gen
         % for 15-64 year old
PrevGeS(zz,ge)=sum(AgeS(zz,4:13,ge))./sum(PopAge(zz,4:13,ge));
PrevGeO(zz,ge)=sum(AgeO(zz,4:13,ge))./sum(PopAge(zz,4:13,ge));
PrevGeF(zz,ge)=sum(AgeF(zz,4:13,ge))./sum(PopAge(zz,4:13,ge));
%for 20-79 year old
PrevGeS2079(zz,ge)=sum(AgeS(zz,5:16,ge))./sum(PopAge(zz,5:16,ge));
PrevGeO2079(zz,ge)=sum(AgeO(zz,5:16,ge))./sum(PopAge(zz,5:16,ge));
PrevGeF2079(zz,ge)=sum(AgeF(zz,5:16,ge))./sum(PopAge(zz,5:16,ge));
    for rfs=1:16
       Prevge(zz,rfs,ge)=sum(PrevAgeGet(zz,rfs,4:13,ge)); 
       
       Prevge2079(zz,rfs,ge)=sum(PrevAgeGet(zz,rfs,5:16,ge));
    end
              
     end
end
%for 15-64 year ald
Prev=sum(Prevge,3);
 %for 20-79 year old
Prev2079=sum(Prevge2079,3);

PrevSALL= (Prev(:,3)+ Prev(:,5)+ Prev(:,7)+Prev(:,8)+Prev(:,11)+ Prev(:,13)+Prev(:,15)+Prev(:,16)); %(S+OS+SF+OSF+SDM+OSDM+SFDM+OSFDM)/NN

PrevOALL=(Prev(:,2)+ Prev(:,5)+ Prev(:,6)+Prev(:,8)+Prev(:,10)+ Prev(:,13)+Prev(:,14)+Prev(:,16)); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN
PrevOALLGE=squeeze(Prevge2079(:,2,:)+ Prevge2079(:,5,:)+ Prevge2079(:,6,:)+Prevge2079(:,8,:)+Prevge2079(:,10,:)+ Prevge2079(:,13,:)+Prevge2079(:,14,:)+Prevge2079(:,16,:)); %(O+OS+OF+OSF+ODM+OSDM+OFDM+OSFDM)/NN

PrevFALL=(Prev(:,4)+ Prev(:,6)+ Prev(:,7)+Prev(:,8)+Prev(:,12)+ Prev(:,14)+Prev(:,15)+Prev(:,16)); % (F+SF+OF+OSF+FDM+OFDM+SFDM+OSFDM)/NN

PrevSum=sum(Prev,2);%S+PrevO+PrevH+PrevOS+PrevOF+PrevSF+PrevOSF
%%DM
%% 2012
%OVERALL
PropDMF12=(sum(DMALLAge(crossdate1,4:13,1)))/(sum(PopAge(crossdate1,4:13,1)));
PropDMM12=(sum(DMALLAge(crossdate1,4:13,2)))/(sum(PopAge(crossdate1,4:13,2)));
% 75+
PropDMFbelo24=(sum(DMALLAge(crossdate1,4:5,1)))/(sum(PopAge(crossdate1,4:5,1)));
PropDMMbelo24=(sum(DMALLAge(crossdate1,4:5,2)))/(sum(PopAge(crossdate1,4:5,2)));

PropDMF25to34=(sum(DMALLAge(crossdate1,6:7,1)))/(sum(PopAge(crossdate1,6:7,1)));
PropDMM25to34=(sum(DMALLAge(crossdate1,6:7,2)))/(sum(PopAge(crossdate1,6:7,2)));

PropDMF35to44=(sum(DMALLAge(crossdate1,8:9,1)))/(sum(PopAge(crossdate1,8:9,1)));
PropDMM35to44=(sum(DMALLAge(crossdate1,8:9,2)))/(sum(PopAge(crossdate1,8:9,2)));

PropDMF45to54=(sum(DMALLAge(crossdate1,10:11,1)))/(sum(PopAge(crossdate1,10:11,1)));
PropDMM45to54=(sum(DMALLAge(crossdate1,10:11,2)))/(sum(PopAge(crossdate1,10:11,2)));

PropDMF55to64=(sum(DMALLAge(crossdate1,12:13,1)))/(sum(PopAge(crossdate1,12:13,1)));
PropDMM55to64=(sum(DMALLAge(crossdate1,12:13,2)))/(sum(PopAge(crossdate1,12:13,2)));

%% -------------------------------------------------------------------------------------- %%
%% RISK FACTORS
%%2012 OBESITY
PropOF12=(sum(AgeO(crossdate1,4:13,1)))/(sum(PopAge(crossdate1,4:13,1)));
PropOM12=(sum(AgeO(crossdate1,4:13,2)))/(sum(PopAge(crossdate1,4:13,2)));
% 75+
PropOFbelo24=(sum(AgeO(crossdate1,4:5,1)))/(sum(PopAge(crossdate1,4:5,1)));
PropOMbelo24=(sum(AgeO(crossdate1,4:5,2)))/(sum(PopAge(crossdate1,4:5,2)));

PropOF25to34=(sum(AgeO(crossdate1,6:7,1)))/(sum(PopAge(crossdate1,6:7,1)));
PropOM25to34=(sum(AgeO(crossdate1,6:7,2)))/(sum(PopAge(crossdate1,6:7,2)));

PropOF35to44=(sum(AgeO(crossdate1,8:9,1)))/(sum(PopAge(crossdate1,8:9,1)));
PropOM35to44=(sum(AgeO(crossdate1,8:9,2)))/(sum(PopAge(crossdate1,8:9,2)));

PropOF45to54=(sum(AgeO(crossdate1,10:11,1)))/(sum(PopAge(crossdate1,10:11,1)));
PropOM45to54=(sum(AgeO(crossdate1,10:11,2)))/(sum(PopAge(crossdate1,10:11,2)));

PropOF55to64=(sum(AgeO(crossdate1,12:13,1)))/(sum(PopAge(crossdate1,12:13,1)));
PropOM55to64=(sum(AgeO(crossdate1,12:13,2)))/(sum(PopAge(crossdate1,12:13,2)));

%%2012 SMOKING
PropSF12=(sum(AgeS(crossdate1,4:13,1)))/(sum(PopAge(crossdate1,4:13,1)));
PropSM12=(sum(AgeS(crossdate1,4:13,2)))/(sum(PopAge(crossdate1,4:13,2)));
% 75+
PropSFbelo24=(sum(AgeS(crossdate1,4:5,1)))/(sum(PopAge(crossdate1,4:5,1)));
PropSMbelo24=(sum(AgeS(crossdate1,4:5,2)))/(sum(PopAge(crossdate1,4:5,2)));

PropSF25to34=(sum(AgeS(crossdate1,6:7,1)))/(sum(PopAge(crossdate1,6:7,1)));
PropSM25to34=(sum(AgeS(crossdate1,6:7,2)))/(sum(PopAge(crossdate1,6:7,2)));

PropSF35to44=(sum(AgeS(crossdate1,8:9,1)))/(sum(PopAge(crossdate1,8:9,1)));
PropSM35to44=(sum(AgeS(crossdate1,8:9,2)))/(sum(PopAge(crossdate1,8:9,2)));

PropSF45to54=(sum(AgeS(crossdate1,10:11,1)))/(sum(PopAge(crossdate1,10:11,1)));
PropSM45to54=(sum(AgeS(crossdate1,10:11,2)))/(sum(PopAge(crossdate1,10:11,2)));

PropSF55to64=(sum(AgeS(crossdate1,12:13,1)))/(sum(PopAge(crossdate1,12:13,1)));
PropSM55to64=(sum(AgeS(crossdate1,12:13,2)))/(sum(PopAge(crossdate1,12:13,2)));

%%2012 SMOKING
PropFF12=(sum(AgeF(crossdate1,4:13,1)))/(sum(PopAge(crossdate1,4:13,1)));
PropFM12=(sum(AgeF(crossdate1,4:13,2)))/(sum(PopAge(crossdate1,4:13,2)));
% 75+
PropFFbelo24=(sum(AgeF(crossdate1,4:5,1)))/(sum(PopAge(crossdate1,4:5,1)));
PropFMbelo24=(sum(AgeF(crossdate1,4:5,2)))/(sum(PopAge(crossdate1,4:5,2)));

PropFF25to34=(sum(AgeF(crossdate1,6:7,1)))/(sum(PopAge(crossdate1,6:7,1)));
PropFM25to34=(sum(AgeF(crossdate1,6:7,2)))/(sum(PopAge(crossdate1,6:7,2)));

PropFF35to44=(sum(AgeF(crossdate1,8:9,1)))/(sum(PopAge(crossdate1,8:9,1)));
PropFM35to44=(sum(AgeF(crossdate1,8:9,2)))/(sum(PopAge(crossdate1,8:9,2)));

PropFF45to54=(sum(AgeF(crossdate1,10:11,1)))/(sum(PopAge(crossdate1,10:11,1)));
PropFM45to54=(sum(AgeF(crossdate1,10:11,2)))/(sum(PopAge(crossdate1,10:11,2)));

PropFF55to64=(sum(AgeF(crossdate1,12:13,1)))/(sum(PopAge(crossdate1,12:13,1)));
PropFM55to64=(sum(AgeF(crossdate1,12:13,2)))/(sum(PopAge(crossdate1,12:13,2)));

%% FITTING OUTPUT DATA
outy=[PropDMFbelo24 PropDMF25to34 PropDMF35to44 PropDMF45to54 PropDMF55to64 PropDMF12...
      PropDMMbelo24 PropDMM25to34 PropDMM35to44 PropDMM45to54 PropDMM55to64 PropDMM12...
      PropOFbelo24 PropOF25to34 PropOF35to44 PropOF45to54 PropOF55to64 PropOF12...
      PropOMbelo24 PropOM25to34 PropOM35to44 PropOM45to54 PropOM55to64 PropOM12...
      PropSFbelo24 PropSF25to34 PropSF35to44 PropSF45to54 PropSF55to64 PropSF12...
      PropSMbelo24 PropSM25to34 PropSM35to44 PropSM45to54 PropSM55to64 PropSM12...
      PropFFbelo24 PropFF25to34 PropFF35to44 PropFF45to54 PropFF55to64 PropFF12...
      PropFMbelo24 PropFM25to34 PropFM35to44 PropFM45to54 PropFM55to64 PropFM12]'.*100; % PHIA
end