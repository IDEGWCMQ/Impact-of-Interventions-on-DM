%% DM MODEL
%% UPDATED January 2021
%% SUSANNE F. AWAD
% clear all
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
Estimates2 = fittingparam2;
Estimates1=fittingparam1;
%%
[X,IncidenceTOT, PARSALL220_79, PARSALL2ge20_79, PAROALL2ge20_79,PAROALL220_79, PARFALL220_79, PARFALL2ge20_79,PrevGeSALL20_79,PrevGeOALL20_79,PrevGeFALL20_79,PrevGeS2079,PrevGeO2079,PrevGeF2079,FALL,OALL,SALL,DMGE2079, PrevDMGE2079, PrevDMALL20_79, S1sumF, S1sumM,PrevDMALL15_64,PrevDMGE,IncidenceRFGE15_64,Nagetot,DMGE,DMNGE,IncidenceRFGE15_64sex,DMAGE2079,IncRate,IncRatesex, PrevOALL, PrevGeO, PrevSALL, PrevGeS, PrevFALL, PrevGeF, PAROALL2, PAROALL2ge, PARSALL2, PARSALL2ge, PARFALL2, PARFALL2ge]=DM_model_3RFs_FINAL(Estimates1,Estimates2, data1, BFf, BMm,Btot);
% % part1F=PAROALLAge2ex(:,:,1);
% % part1F1=PAROALLAge2ex(:,:,2);
close all
%% FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
options=optimset('Display','iter', 'MaxFunEvals',1000,'MaxIter',1000, 'TolFun', 1e-5, 'TolX',1e-4);
%% ESTIMATED PREVALENCES 
% [outy]=DM_model_3RFsDM(fittingparam1,Result);
% %% upper and lower bounds
% %females
% lb(1)=0.00000001;
% lb(2)=0.1550;
% lb(3)=0.01;
% 
% ub(1)=1;
% ub(2)=0.2100;
% ub(3)=10000;
% 
% 
% lb(4)=0.00000001;
% lb(5)=0.00000001;
% 
% ub(4)=2.0;
% ub(5)=10000;
% 
% 
% lb(6)=0.00000001;
% lb(7)=0.1550;
% lb(8)=0.01;
% 
% ub(6)=1;
% ub(7)=0.2500;
% ub(8)=10000;
% %MAles 
% lb(9)=0.00000001;
% lb(10)=0.1550;
% lb(11)=0.01;
% 
% ub(9)=1;
% ub(10)=0.2100;
% ub(11)=10000;
% 
% 
% lb(12)=0.00000001;
% lb(13)=0.00000001;
% 
% ub(12)=2.0;
% ub(13)=10000;
% 
% 
% lb(14)=0.00000001;
% lb(15)=0.1550;
% lb(16)=0.01;
% 
% ub(14)=1;
% ub(15)=0.2500;
% ub(16)=10000;
% %females
% lb(17)=[0];
% lb(18)=0;
% % lb(19)=0.19;
% 
% ub(17)=20;
% ub(18)=20;
% % ub(19)=0.4;

%% Fiting DM prevalence
%fittingparam1=Estimates1;
lb=ones(120,1)*10^(-100000);
ub=ones(120,1).*2;

lb(121:132,1)=1*10^(-1000);
ub(121:132,1)=1*500;

[Final_param1,fval,exitflag,output] = fminsearchbnd(@DMminfun,fittingparam1,lb,ub,options,Result,data1, BFf, BMm,Btot);
Estimates1=Final_param1;% Estimates1=fittingparam1
% Final_param1= fittingparam1;
% % Fitting population
options=optimset('Display','iter', 'MaxFunEvals',2000,'MaxIter',2000, 'TolFun', 1e-5, 'TolX',1e-4);
 for xixxx=1:2
     fittingparam2=Estimates2;
[Final_param2,fval,exitflag,output] = fminsearchbnd(@minfun,fittingparam2,[],[],options,Final_param1,data1, BFf, BMm,Btot);%[],[]lbub
Estimates2 = Final_param2;
 end
% %  fittingparam2=Estimates2;
toc
%% ESTIMATED PREVALENCES 
close all
[S1sumF, S1sumM,alphasn, p1, p2, PrevDMALL15_64,PrevDMGE,IncidenceRFGE15_64,Nagetot,cost_per_case2f,cost_per_case3f,cost_per_case2l,cost_per_case3l,DMGE,DMNGE,IncidenceRFGE15_64sex,DMAGE2079,IncRate,IncRatesex, PrevOALL, PrevGeO, PrevSALL, PrevGeS, PrevFALL, PrevGeF, PAROALL2, PAROALL2ge, PARSALL2, PARSALL2ge, PARFALL2, PARFALL2ge]=DM_model_3RFs_FINAL(Estimates1,Estimates2, data1, BFf, BMm,Btot);
% [S1sumF, S1sumM]=DM_model_3RFs_FINALDEMO(Estimates1,Estimates2);

%% Redefine parameter

finalE=Estimates1(1:20*2*3,1);
param=reshape(finalE,20,2,3);
paramDM=Estimates1(20*2*3+1:(20*2*3)+(6),1);
paramOB=Estimates1((20*2*3)+(6)+1:(20*2*3+6)+(3*2),1);
% Anch=Estimates1(139:141,1);
ResultAll=Estimates2;%(1:8,1);
% Result(1:8,2)=Estimates2(9:16,1);
% Result(9,1)=Estimates2(17,1);
% Result(9,2)=Estimates2(18,1);
% Result(10,1)=Estimates2(19,1);
% Result(10,2)=Estimates2(20,1);% Result(12:14,1)=Estimates2(23:25,1);
% Result(12:14,2)=Estimates2(26:28,1);


% Result(1:8,1)=Estimates2(1:8,1);
% Result(1:8,2)=Estimates2(9:16,1);

save('pram1DEC17.mat', 'param')
save('pramDMDEC17.mat', 'paramDM')
save('pramOBDEC17.mat', 'paramOB')
save('ResultQATARnew.mat', 'ResultAll')
% save('Anchor.mat','Anch')
save('parameterRates.mat', 'Estimates1')
save('parameterPop.mat', 'Estimates2')
% save('prevalence.mat', 'f_prev')
% %%
save