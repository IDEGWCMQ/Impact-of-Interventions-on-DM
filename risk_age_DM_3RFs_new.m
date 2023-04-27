function [dxL,alpha,dmmort]=risk_age_DM_3RFs_new(t,xt,lambda,eta,muALL,RRm,alphas, beta, shei,sigma,scenario,Result,lambdaf,lambdafi,lambdaH,lambdaO,lambdaS,lambdaF,lambdaOS,lambdaOF,lambdaSF, lambdaOSF,intv)

global ag gen NcompTOT N0 rf NrfTOT delta phi epsilon iota q tt d a rho theta pii ii p ss b cc nu xi sampi w n hh e etaa omega kappa xx ff gg oo gamma upsilon r u k iii Leps chi psi Omg m jj v l 
alpha=ones(ag,gen);
alphapart1=(OBincFuncGaussian(ag,gen,alphas(1,:),alphas(2,:),alphas(3,:)));
% if t<=2017
for gender=1:gen
    for aging=1:ag
    alpha(aging,gender)=alphapart1(aging,gender);%*(alphas(4,gender)./(1 + exp(-alphas(5,gender)*(t-alphas(6,gender)))));
    end
end
alpha(1:3,:)=0;
alpha(4,:)=0.005;
%becoming obese
chi=alpha;%0;%
psi=alpha;%0;%
Omg=alpha;%0;%
m=alpha;%alpha;
jj=alpha;%alpha;
v=alpha;%alpha;
l=alpha;%alpha;
% scenario=3.1;
%%
% dmmort=zeros(1,ag,gen);

x=reshape(xt,NcompTOT*intv,ag,gen);
N=squeeze(sum(sum(x,2),1));

%% LSM
IEH=ones(8,20);
da=zeros(NcompTOT*intv,ag,gen);
tv=2056;
scale2=0;
scale3=0;
scale1=0;
if scenario==1
    IEH=ones(8,20);
    da=zeros(NcompTOT*intv,ag,gen);
tv=2056;
scale2=0;
scale3=0;
scale1=0;
end
if scenario>3 && scenario<4
    if scenario==3.1 || scenario==3.6; IEH(2,:)=0.70; IEH(5,:)=0.70; IEH(6,:)=0.70; IEH(8,:)=0.70; % OBESITY + ANY RF (REGARDLESS OF AGE)-50%COVERAGE LSM EFFECACY=0.65
% elseif scenario==3.2 || scenario==3.7; IEH(2,:)=0.79; IEH(5,:)=0.79; IEH(6,:)=0.79; IEH(8,:)=0.79; % OBESITY + ANY RF (REGARDLESS OF AGE)-50%COVERAGE LSM EFFECACY=0.79
elseif scenario==3.3 || scenario==3.8; IEH(2,8:20)=0.70; IEH(5,8:20)=0.70; IEH(6,8:20)=0.70; IEH(8,8:20)=0.70;% OBESITY + AGE 35+ -50%COVERAGE LSM EFFECACY=0.65
% elseif scenario==3.4 || scenario==3.9; IEH(2,8:20)=0.79; IEH(5,8:20)=0.79; IEH(6,8:20)=0.79; IEH(8,8:20)=0.79;% OBESITY + AGE 35+ -50%COVERAGE LSM EFFECACY=0.79
elseif scenario==3.5 || scenario==3.11; IEH(1:8,11:20)=0.70;
elseif scenario==3.12 || scenario==3.13; IEH(1:8,11:20)=0.60;

end      
lambdafi=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];
lambda(:,:,1)=[lambdaf(1:2:end,:)
    zeros(8,20)
    lambdafi(1:2:end,:)
     zeros(8,20)] ;%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=[lambdaf(2:2:end,:)
    zeros(8,20)
    lambdafi(2:2:end,:)
    zeros(8,20)];%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,1,:)=0;

if t>2021 && t<=2026
    if scenario==3.1; scale1=0.27; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(8,:,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.2; scale1=0.27; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(8,:,:)=1.*(scale1+scale2+scale3); 
    elseif scenario==3.3; scale1=0.95; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.4; scale1=0.95; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.5; scale1=1000; scale2=0; scale3=0; da(1:8,11:20,:)=1.*scale1; %da(2,11:20,:)=1.*scale1; da(5,11:20,:)=1.*(scale1+scale2); da(6,11:20,:)=1.*(scale1+scale3); da(3,11:20,:)=1.*scale1+scale2; da(4,11:20,:)=1.*scale1+scale3; da(7,11:20,:)=1.*(scale1+scale2+scale3); da(8,11:20,:)=1.*(scale1+scale2+scale3);

    elseif scenario==3.6; scale1=0.109; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,:,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.7; scale1=0.109; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,:,:)=1.*(scale1+scale2+scale3); 
    elseif scenario==3.8; scale1=0.205; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.9; scale1=0.205; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
    
    elseif scenario==3.11;scale1=0.75; scale2=0; scale3=0; da(1:8,11:20,:)=1.*scale1; %da(2,11:20,:)=1.*scale1; da(5,11:20,:)=1.*(scale1+scale2); da(6,11:20,:)=1.*(scale1+scale3); da(3,11:20,:)=1.*scale1+scale2; da(4,11:20,:)=1.*scale1+scale3; da(7,11:20,:)=1.*(scale1+scale2+scale3); da(8,11:20,:)=1.*(scale1+scale2+scale3);
    
%     else; da=zeros(NcompTOT*intv,ag,gen);
    end
elseif t>2026
    if scenario==3.1; scale1=0.009; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,:,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.2; scale1=0.009; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,:,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.3; scale1=0.019; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3); 
    elseif scenario==3.4; scale1=0.019; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);

    elseif scenario==3.6; scale1=0.003; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,:,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.7; scale1=0.003; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,:,:)=1.*(scale1+scale2+scale3); 
    elseif scenario==3.8; scale1=0.005; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.9; scale1=0.005; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
    
    elseif scenario==3.5; scale1=0.15; scale2=0; scale3=0; da(1:8,11:20,:)=1.*scale1; %da(2,11:20,:)=1.*scale1; da(5,11:20,:)=1.*(scale1+scale2); da(6,11:20,:)=1.*(scale1+scale3); da(3,11:20,:)=1.*scale1+scale2; da(4,11:20,:)=1.*scale1+scale3; da(7,11:20,:)=1.*(scale1+scale2+scale3); da(8,11:20,:)=1.*(scale1+scale2+scale3);
    elseif scenario==3.11;scale1=0.019; scale2=0; scale3=0; da(1,11:20,:)=1.*scale1; da(2,11:20,:)=1.*scale1; da(5,11:20,:)=1.*(scale1+scale2); da(6,11:20,:)=1.*(scale1+scale3); da(3,11:20,:)=1.*scale1+scale2; da(4,11:20,:)=1.*scale1+scale3; da(7,11:20,:)=1.*(scale1+scale2+scale3); da(8,11:20,:)=1.*(scale1+scale2+scale3);

 
    end
    
end
end
%     else; %da=zeros(NcompTOT*intv,ag,gen);
% end
% TRANSPORTATION
if scenario>4 && scenario<5
if scenario==4.1 || scenario==4.2; IEH(1:8,:)=0.75; % 
elseif scenario==4.3 || scenario==4.4  || scenario==4.5 || scenario==4.6; IEH(4,:)=1/1.48; IEH(7,:)=1/1.48; IEH(6,:)=1/1.48; IEH(8,:)=1/1.48; % OBESITY + ANY RF (REGARDLESS OF AGE)-50%COVERAGE LSM EFFECACY=0.79  
end
    lambdafi=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];

lambda(:,:,1)=[lambdaf(1:2:end,:)
    zeros(8,20)
    lambdafi(1:2:end,:)
     zeros(8,20)] ;%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=[lambdaf(2:2:end,:)
    zeros(8,20)
    lambdafi(2:2:end,:)
    zeros(8,20)];%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,1,:)=0;

if scenario==4.1 || scenario==4.3 || scenario==4.5; tv=2030;
elseif scenario==4.2 || scenario==4.4 || scenario==4.6; tv=2031;
end
if t>2021 && t<=tv
    if scenario==4.1; scale1=0.0113; elseif scenario==4.2; scale1=0.011; elseif scenario==4.3 
    scale1=0.0097; elseif scenario==4.4; scale1=0.012; elseif scenario==4.5; scale1=0.0125; 
    elseif scenario==4.6; scale1=0.013; end
da(1:8,:,:)=1.*scale1; 
elseif t>tv
if scenario==4.1; scale1=0.0032; elseif scenario==4.2; scale1=0.003; elseif scenario==4.3 
   scale1=0.0032; elseif scenario==4.4; scale1=0.0033; elseif scenario==4.5; scale1=0.0034; 
    elseif scenario==4.6; scale1=0.0033; end
da(1:8,:,:)=1.*scale1;
% else
%   da(2:8,:,:) =0;   
end
 end
%% LEGISLATION
if scenario>5 && scenario<6
if scenario==5.1 || scenario==5.3 || scenario==5.5 || scenario==5.6; IEH(1:8,:)=1 ; elseif scenario==5.2 || scenario==5.4 || scenario==5.7 || scenario==5.8; IEH(1:8,:)=0.93; % 
end
    lambdafi=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];

lambda(:,:,1)=[lambdaf(1:2:end,:)
    zeros(8,20)
    lambdafi(1:2:end,:)
     zeros(8,20)] ;%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=[lambdaf(2:2:end,:)
    zeros(8,20)
    lambdafi(2:2:end,:)
    zeros(8,20)];%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,1,:)=0;

if scenario==5.1 || scenario==5.2 || scenario==5.5 || scenario==5.7; tv=2024;
elseif scenario==5.3 || scenario==5.4 || scenario==5.6 || scenario==5.8; tv=2050;
end
if t>2021 && t<=tv
    if scenario==5.1; scale1=0.0105; elseif scenario==5.2; scale1=0.0105; elseif scenario==5.3 
    scale1=0.006; elseif scenario==5.4; scale1=0.006; elseif scenario==5.5; scale1=0.0125; 
    elseif scenario==5.6; scale1=0.007; elseif scenario==5.7; scale1=0.021; elseif scenario==5.8; scale1=0.009;
    end
% else                % PA:%0.009 no change;
da(1:8,:,:)=1.*scale1;
elseif t>tv
if scenario==5.1; scale1=0.0032; elseif scenario==5.2; scale1=0.0032; elseif scenario==5.3 
   scale1=0.0035; elseif scenario==5.4; scale1=0.0035; elseif scenario==5.5; scale1=0.0034; 
    elseif scenario==5.6; scale1=0.0033; elseif scenario==5.7; scale1=0.0035; 
end
da(1:8,:,:)=1.*scale1;   
end
end
%% CONSUMPSION
if scenario>6 && scenario<7
    if scenario==6.1 || scenario==6.6; IEH(1:8,:)=0.93; lambdafi=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];
elseif scenario==6.2 || scenario==6.7; IEH(1:8,:)=0.90; lambdafi=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];
elseif scenario==6.3 || scenario==6.7; IEH(1:8,:)=0.87; lambdafi=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];
elseif scenario==6.4 || scenario==6.8; IEH(1:8,:)=1;    lambdafi=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
elseif scenario==6.5 || scenario==6.9; IEH(1:8,:)=0.93; lambdafi=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
end      
lambda(:,:,1)=[lambdaf(1:2:end,:)
    zeros(8,20)
    lambdafi(1:2:end,:)
     zeros(8,20)] ;%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=[lambdaf(2:2:end,:)
    zeros(8,20)
    lambdafi(2:2:end,:)
    zeros(8,20)];%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,1,:)=0;

if scenario==6.1 || scenario==6.2 || scenario==6.3 || scenario==6.7; tv=2026;
elseif scenario==6.5 || scenario==6.4 || scenario==6.6 || scenario==6.8; tv=2024;
end

if t>2021 && t<=tv
    if scenario==6.1 || scenario==6.2 || scenario==6.3; scale1=0.22; scale2=0; scale3=0; da(:,:,:)=1.*scale1; 
elseif scenario==6.4; scale1=0.12; scale2=0; scale3=0; da(:,4:13,:)=1.*scale1;
elseif scenario==6.5; scale1=0.12; scale2=0; scale3=0; da(:,4:13,:)=1.*scale1;
elseif scenario==6.6; scale1=0.109; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,:,:)=1.*(scale1+scale2+scale3); 
elseif scenario==6.7; scale1=0.205; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
elseif scenario==6.8; scale1=0.205; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
    
elseif scenario==6.9; scale1=1000; scale2=0; scale3=0; da(1:8,11:20,:)=1.*scale1; %da(2,11:20,:)=1.*scale1; da(5,11:20,:)=1.*(scale1+scale2); da(6,11:20,:)=1.*(scale1+scale3); da(3,11:20,:)=1.*scale1+scale2; da(4,11:20,:)=1.*scale1+scale3; da(7,11:20,:)=1.*(scale1+scale2+scale3); da(8,11:20,:)=1.*(scale1+scale2+scale3);
elseif scenario==6.11;scale1=0.75; scale2=0; scale3=0; da(1:8,11:20,:)=1.*scale1; %da(2,11:20,:)=1.*scale1; da(5,11:20,:)=1.*(scale1+scale2); da(6,11:20,:)=1.*(scale1+scale3); da(3,11:20,:)=1.*scale1+scale2; da(4,11:20,:)=1.*scale1+scale3; da(7,11:20,:)=1.*(scale1+scale2+scale3); da(8,11:20,:)=1.*(scale1+scale2+scale3);
    end
    
elseif t>tv
    if scenario==6.1 || scenario==6.2 || scenario==6.3; scale1=0.03; scale2=0; scale3=0; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1+scale2); da(6,:,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,:,:)=1.*(scale1+scale2+scale3);
elseif scenario==6.4; scale1=0.015; scale2=0; scale3=0; da(:,4:13,:)=1.*scale1;
elseif scenario==6.5; scale1=0.008105; scale2=0; scale3=0; da(:,4:13,:)=1.*scale1;
    end
end
end% da(1:8,:,:)=1.*scale1;
% %% combination
if scenario==7
    if scenario==7; IEH(2,8:20)=0.75*0.47; IEH(5,8:20)=0.75*0.47; IEH(6,8:20)=0.75*0.47; IEH(8,8:20)=0.75*0.47;
        IEH(1:8,4:13)=0.75*0.65; IEH(1,:)=0.75*0.70;IEH(3,:)=0.75*0.70;IEH(4,:)=0.75*0.70; IEH(7,:)=0.75*0.70;% OBESITY + AGE 35+ -50%COVERAGE LSM EFFECACY=0.65
end      
    lambdafi=[lambdaH.*IEH(1,:); lambdaH.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaS.*IEH(5,:); lambdaF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaSF.*IEH(8,:)];
lambda(:,:,1)=[lambdaf(1:2:end,:)
    zeros(8,20)
    lambdafi(1:2:end,:)
     zeros(8,20)] ;%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=[lambdaf(2:2:end,:)
    zeros(8,20)
    lambdafi(2:2:end,:)
    zeros(8,20)];%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,1,:)=0;

if t>2021 && t<=2026
%     if scenario==7 
    scale1N=0.019; da(:,1:7,:)=1.*scale1N; da(1,:,:)=1.*scale1N; da(3,:,:)=1.*scale1N; da(4,:,:)=1.*scale1N; da(7,:,:)=1.*scale1N;
    scale1=0.95; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
%     end
elseif t>2026
%     if scenario==7;
        scale1=0.019; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3); 
    scale1N=0.009; da(:,1:7,:)=1.*scale1N; da(1,:,:)=1.*scale1N; da(3,:,:)=1.*scale1N; da(4,:,:)=1.*scale1N; da(7,:,:)=1.*scale1N;

%     else; da=zeros(NcompTOT*intv,ag,gen); 
%     end
end
end
if scenario==8
    IEH(1:8,:)=.93; IEH(2,:)=0.70*.93; IEH(5,:)=0.70*.93; IEH(6,:)=0.70*.93; IEH(8,:)=0.70*.93; 
    lambdafi=[lambdaH.*IEH(1,:); lambdaO.*IEH(2,:); lambdaS.*IEH(3,:); lambdaF.*IEH(4,:); lambdaOS.*IEH(5,:); lambdaOF.*IEH(6,:); lambdaSF.*IEH(7,:); lambdaOSF.*IEH(8,:)];
lambda(:,:,1)=[lambdaf(1:2:end,:)
    zeros(8,20)
    lambdafi(1:2:end,:)
     zeros(8,20)] ;%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,:,2)=[lambdaf(2:2:end,:)
    zeros(8,20)
    lambdafi(2:2:end,:)
    zeros(8,20)];%reshape(lambdaf,NrfTOT+1,length(lambdaf),gen);
lambda(:,1,:)=0;

if t>2021 && t<=2026
    scale1=0.097; scale2=0.01; scale3=0; 
    da(1:8,:,:)=1.*scale2; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1); da(6,:,:)=1.*(scale1); da(8,:,:)=1.*(scale1+scale2+scale3);
    
%     scale1N=0.009; da(:,1:7,:)=1.*scale1N; da(1,:,:)=1.*scale1N; da(3,:,:)=1.*scale1N; da(4,:,:)=1.*scale1N; da(7,:,:)=1.*scale1N;
%     scale1=0.95; scale2=0; scale3=0; da(2,8:20,:)=1.*scale1; da(5,8:20,:)=1.*(scale1+scale2); da(6,8:20,:)=1.*(scale1+scale3); da(3,:,:)=1.*scale2; da(4,:,:)=1.*scale3; da(7,:,:)=1.*(scale2+scale3); da(8,8:20,:)=1.*(scale1+scale2+scale3);
%     end
elseif t>2026
    scale1=0.008; scale2=0.01; scale3=0; 
    da(1:8,:,:)=1.*scale2; da(2,:,:)=1.*scale1; da(5,:,:)=1.*(scale1); da(6,:,:)=1.*(scale1); da(8,:,:)=1.*(scale1+scale2+scale3);

%     else; da=zeros(NcompTOT*intv,ag,gen); 
    end
 end
% end
%%
% Gaussian distribution Birth for t,2005

a1f=(Result(1,1));%[0.0910590437029461];%squeeze[0.0376398585222292];%0.04103620599;%
b1f=(Result(2,1))*10^4;%[0.169496261886811]*10^4;%squeeze[0.191619149119278]*10^4;%0.193046793263662*10^4;%;
c1f=(Result(3,1));%[316.534154487845];%squeeze[125.838356613016];%121.07936889569;%;

a1m=(Result(1,2));%[0.105956927135559];%squeeze%[0.0376398585222292];%0.04103620599;%
b1m=(Result(2,2))*10^4;%[0.169761224452025]*10^4;%squeeze[0.191619149119278]*10^4;%0.193046793263662*10^4;%;
c1m=(Result(3,2));%[273.094561242040];%squeeze[125.838356613016];%121.07936889569;%;
% % Mortality rate
% % the gaussian part
b2f=Result(4,:).*10;%[1.38823493843683]*10;%14;%
c2f=Result(5,:);%[1.60106043463158];%[1.64822373335293];%
% the logistic part
a3f=Result(6,:);%[0.129148406527900];%0.14300;%[0.119852332457220];%
b3f=Result(7,:).*10^4;%[0.182453620234226]*10^4;%[0.200418906213193]*10^4;%
c3f=Result(8,:);

s=Result(9,:);% b4f=y(10)*10^4;
%%
dx=zeros(NcompTOT*intv,ag,gen);
% dmmort=zeros(1,20,2);
if t<1930
%     Immigration(1:ag,1:gen)=0;
dx(1,1,1) =N(1,1).*BirthFuncGaussian(t,a1f,b1f,c1f)-Mortality(t,1,b2f(1,1),c2f(1,1),a3f(1,1),b3f(1,1),c3f(1,1),s(1,1))*x(1,1,1)-eta(1)*x(1,1,1); 
dx(1,1,2) =N(2,1).*BirthFuncGaussian(t,a1m,b1m,c1m)-Mortality(t,1,b2f(1,2),c2f(1,2),a3f(1,2),b3f(1,2),c3f(1,2),s(1,2))*x(1,1,2)-eta(1)*x(1,1,2); 
% dmmort(1,1)=0;
% dmmort(1,2)=0;
% mmort(1,1)=Mortality(t,1,b2f(1,1),c2f(1,1),a3f(1,1),b3f(1,1),c3f(1,1),s(1,1))*x(1,1,1);
% mmort(1,2)=Mortality(t,1,b2f(1,2),c2f(1,2),a3f(1,2),b3f(1,2),c3f(1,2),s(1,2))*x(1,1,2);

for agn=2:ag
    for ge=1:gen
dx(1,agn,ge) =-Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(1,agn,ge)+eta(agn-1)*x(1,agn-1,ge)-eta(agn)*x(1,agn,ge); %*brate*N*mu(agn,ge)*N0
% dmmort(agn,ge)=0;
for gr=1:32
dmmort(1,agn,ge,gr)=Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(gr,agn,ge);
end
    end
end
elseif t>=2005 && t<2023
    
dx(1,1,1) =N(1,1).*BirthFuncGaussian(t,a1f,b1f,c1f)-Mortality(t,1,b2f(1,1),c2f(1,1),a3f(1,1),b3f(1,1),c3f(1,1),s(1,1))*x(1,1,1)-eta(1)*x(1,1,1); %*brate*N*mu(agn,ge)*N0
dx(1,1,2) =N(2,1).*BirthFuncGaussian(t,a1m,b1m,c1m)-Mortality(t,1,b2f(1,2),c2f(1,2),a3f(1,2),b3f(1,2),c3f(1,2),s(1,2))*x(1,1,2)-eta(1)*x(1,1,2); %*brate*N*mu(agn,ge)*N0   
% dmmort(1,1)=0;
% dmmort(1,2)=0;
% mmort(1,1)=Mortality(t,1,b2f(1,1),c2f(1,1),a3f(1,1),b3f(1,1),c3f(1,1),s(1,1))*x(1,1,1);
% mmort(1,2)=Mortality(t,1,b2f(1,2),c2f(1,2),a3f(1,2),b3f(1,2),c3f(1,2),s(1,2))*x(1,1,2);

    for agn=2:ag
 for ge=1:gen
%      Immigration(agn,ge)=(a4m(1,ge)*exp(-((agn-b4m(1,ge))/c4m(1,ge))^2))*exp(-((t-b5m(1,ge))/c5m(1,ge))^2);
dx(1,agn,ge) = -da(1,agn,ge)*x(1,agn,ge)+sigma(agn,ge)*x(2,agn,ge)+delta(agn,ge)*x(3,agn,ge)+phi(agn,ge)*x(4,agn,ge)-lambda(1,agn,ge)*x(1,agn,ge)- beta(agn,ge)*x(1,agn,ge)-alpha(agn,ge)*x(1,agn,ge)-shei(agn,ge)*x(1,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(1,agn,ge)+eta(agn-1)*x(1,agn-1,ge)-eta(agn)*x(1,agn,ge);
% {Obese}
dx(2,agn,ge) = -da(2,agn,ge)*x(2,agn,ge)+alpha(agn,ge)*x(1,agn,ge)+ epsilon(agn,ge)*x(5,agn,ge)+theta(agn,ge)*x(6,agn,ge)-lambda(2,agn,ge)*x(2,agn,ge)-nu(agn,ge)*x(2,agn,ge)- etaa(agn,ge)*x(2,agn,ge)-sigma(agn,ge)*x(2,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(2,agn,ge)+eta(agn-1)*x(2,agn-1,ge)-eta(agn)*x(2,agn,ge);
% {Smoker}
dx(3,agn,ge) = -da(3,agn,ge)*x(3,agn,ge)+beta(agn,ge)*x(1,agn,ge)+gamma(agn,ge)*x(5,agn,ge)+pii(agn,ge)*x(7,agn,ge)-lambda(3,agn,ge)*x(3,agn,ge)- chi(agn,ge)*x(3,agn,ge) - delta(agn,ge)*x(3,agn,ge)-omega(agn,ge)*x(3,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(3,agn,ge)+eta(agn-1)*x(3,agn-1,ge)-eta(agn)*x(3,agn,ge);
% {Physical Active};
dx(4,agn,ge) =  -da(4,agn,ge)*x(4,agn,ge)+shei(agn,ge)*x(1,agn,ge)+rho(agn,ge)*x(7,agn,ge)+Leps(agn,ge)*x(6,agn,ge)-lambda(4,agn,ge)*x(4,agn,ge)  - phi(agn,ge)*x(4,agn,ge) - psi(agn,ge)*x(4,agn,ge)-xi(agn,ge)*x(4,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(4,agn,ge)+eta(agn-1)*x(4,agn-1,ge)-eta(agn)*x(4,agn,ge);

% {overlap compartment}
dx(5,agn,ge) = - da(5,agn,ge)*x(5,agn,ge)+nu(agn,ge)*x(2,agn,ge) +chi(agn,ge)*x(3,agn,ge)+ii(agn,ge)*x(8,agn,ge)- lambda(5,agn,ge)*x(5,agn,ge)- epsilon(agn,ge)*x(5,agn,ge)- gamma(agn,ge)*x(5,agn,ge)-kappa(agn,ge)*x(5,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(5,agn,ge)+eta(agn-1)*x(5,agn-1,ge)-eta(agn)*x(5,agn,ge);

dx(6,agn,ge) = - da(6,agn,ge)*x(6,agn,ge)+etaa(agn,ge)*x(2,agn,ge) +psi(agn,ge)*x(4,agn,ge)+iota(agn,ge)*x(8,agn,ge)- lambda(6,agn,ge)*x(6,agn,ge)-Leps(agn,ge)*x(6,agn,ge)-theta(agn,ge)*x(6,agn,ge)-sampi(agn,ge)*x(6,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(6,agn,ge)+eta(agn-1)*x(6,agn-1,ge)-eta(agn)*x(6,agn,ge);

dx(7,agn,ge) =- da(7,agn,ge)*x(7,agn,ge)+ omega(agn,ge)*x(3,agn,ge) +xi(agn,ge)*x(4,agn,ge)+upsilon(agn,ge)*x(8,agn,ge)- lambda(7,agn,ge)*x(7,agn,ge)- pii(agn,ge)*x(7,agn,ge)- rho(agn,ge)*x(7,agn,ge)-Omg(agn,ge)*x(7,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(7,agn,ge)+eta(agn-1)*x(7,agn-1,ge)-eta(agn)*x(7,agn,ge);

dx(8,agn,ge) =- da(8,agn,ge)*x(8,agn,ge)+ kappa(agn,ge)*x(5,agn,ge)+sampi(agn,ge)*x(6,agn,ge)+Omg(agn,ge)*x(7,agn,ge)- lambda(8,agn,ge)*x(8,agn,ge)- iota(agn,ge)*x(8,agn,ge)- ii(agn,ge)*x(8,agn,ge)-upsilon(agn,ge)*x(8,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(8,agn,ge)+eta(agn-1)*x(8,agn-1,ge)-eta(agn)*x(8,agn,ge);

% {Infected DM}
dx(9,agn,ge) = lambda(1,agn,ge)*x(1,agn,ge) + r(agn,ge)*x(10,agn,ge) + q(agn,ge)*x(11,agn,ge)+p(agn,ge)*x(12,agn,ge)- m(agn,ge)*x(9,agn,ge) - n(agn,ge)*x(9,agn,ge)-oo(agn,ge)*x(9,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)+eta(agn-1)*x(9,agn-1,ge)...
    -eta(agn)*x(9,agn,ge);

dx(10,agn,ge) = lambda(2,agn,ge)*x(2,agn,ge) + a(agn,ge)*x(13,agn,ge) + b(agn,ge)*x(14,agn,ge)+m(agn,ge)*x(9,agn,ge)- e(agn,ge)*x(10,agn,ge) - ff(agn,ge)*x(10,agn,ge)-r(agn,ge)*x(10,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)+eta(agn-1)*x(10,agn-1,ge)...
    -eta(agn)*x(10,agn,ge);

dx(11,agn,ge)= lambda(3,agn,ge)*x(3,agn,ge) + n(agn,ge)*x(9,agn,ge)+iii(agn,ge)*x(13,agn,ge)+cc(agn,ge)*x(15,agn,ge) - jj(agn,ge)*x(11,agn,ge) - gg(agn,ge)*x(11,agn,ge) -q(agn,ge)*x(11,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge) -(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge)+eta(agn-1)*x(11,agn-1,ge)...
    -eta(agn)*x(11,agn,ge);

dx(12,agn,ge)= lambda(4,agn,ge)*x(4,agn,ge) + oo(agn,ge)*x(9,agn,ge)+k(agn,ge)*x(14,agn,ge)+ d(agn,ge)*x(15,agn,ge) - p(agn,ge)*x(12,agn,ge) - l(agn,ge)*x(12,agn,ge) - hh(agn,ge)*x(12,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge) -(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge)+eta(agn-1)*x(12,agn-1,ge)...
    -eta(agn)*x(12,agn,ge);

dx(13,agn,ge)= lambda(5,agn,ge)*x(5,agn,ge) + ss(agn,ge)*x(16,agn,ge)+ jj(agn,ge)*x(11,agn,ge) + e(agn,ge)*x(10,agn,ge) - xx(agn,ge)*x(13,agn,ge) - a(agn,ge)*x(13,agn,ge)-iii(agn,ge)*x(13,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)+eta(agn-1)*x(13,agn-1,ge)...
    -eta(agn)*x(13,agn,ge);

dx(14,agn,ge)= lambda(6,agn,ge)*x(6,agn,ge) + tt(agn,ge)*x(16,agn,ge)+l(agn,ge)*x(12,agn,ge) + ff(agn,ge)*x(10,agn,ge)- w(agn,ge)*x(14,agn,ge)-b(agn,ge)*x(14,agn,ge)-k(agn,ge)*x(14,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)+eta(agn-1)*x(14,agn-1,ge)...
    -eta(agn)*x(14,agn,ge);

dx(15,agn,ge)= lambda(7,agn,ge)*x(7,agn,ge) + u(agn,ge)*x(16,agn,ge)+hh(agn,ge)*x(12,agn,ge) + gg(agn,ge)*x(11,agn,ge)- d(agn,ge)*x(15,agn,ge)-cc(agn,ge)*x(15,agn,ge)-v(agn,ge)*x(15,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)+eta(agn-1)*x(15,agn-1,ge)...
    -eta(agn)*x(15,agn,ge);

dx(16,agn,ge)= lambda(8,agn,ge)*x(8,agn,ge) +xx(agn,ge)*x(13,agn,ge) + w(agn,ge)*x(14,agn,ge)+v(agn,ge)*x(15,agn,ge) - ss(agn,ge)*x(16,agn,ge) -  tt(agn,ge)*x(16,agn,ge)- u(agn,ge)*x(16,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge)+eta(agn-1)*x(16,agn-1,ge)...
    -eta(agn)*x(16,agn,ge);

%% invtervention
dx(17,agn,ge) = da(1,agn,ge)*x(1,agn,ge)+sigma(agn,ge)*x(18,agn,ge)+delta(agn,ge)*x(19,agn,ge)+phi(agn,ge)*x(20,agn,ge)-lambda(17,agn,ge)*x(17,agn,ge)- beta(agn,ge)*x(17,agn,ge)-alpha(agn,ge)*x(17,agn,ge)-shei(agn,ge)*x(17,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(17,agn,ge)+eta(agn-1)*x(17,agn-1,ge)-eta(agn)*x(17,agn,ge);
% {Obese}
dx(18,agn,ge) = da(2,agn,ge)*x(2,agn,ge)+alpha(agn,ge)*x(17,agn,ge)+ epsilon(agn,ge)*x(21,agn,ge)+theta(agn,ge)*x(22,agn,ge)-lambda(18,agn,ge)*x(18,agn,ge)-nu(agn,ge)*x(18,agn,ge)-etaa(agn,ge)*x(18,agn,ge)-sigma(agn,ge)*x(18,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(18,agn,ge)+eta(agn-1)*x(18,agn-1,ge)-eta(agn)*x(18,agn,ge);
% {Smoker}
dx(19,agn,ge) = da(3,agn,ge)*x(3,agn,ge)+beta(agn,ge)*x(17,agn,ge)+gamma(agn,ge)*x(21,agn,ge)+pii(agn,ge)*x(23,agn,ge)-lambda(19,agn,ge)*x(19,agn,ge)- chi(agn,ge)*x(19,agn,ge)-delta(agn,ge)*x(19,agn,ge)-omega(agn,ge)*x(19,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(19,agn,ge)+eta(agn-1)*x(19,agn-1,ge)-eta(agn)*x(19,agn,ge);
% {Physical Active};
dx(20,agn,ge) = da(4,agn,ge)*x(4,agn,ge)+shei(agn,ge)*x(17,agn,ge)+rho(agn,ge)*x(23,agn,ge)+Leps(agn,ge)*x(22,agn,ge)-lambda(20,agn,ge)*x(20,agn,ge)- phi(agn,ge)*x(20,agn,ge)-psi(agn,ge)*x(20,agn,ge)-xi(agn,ge)*x(20,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(20,agn,ge)+eta(agn-1)*x(20,agn-1,ge)-eta(agn)*x(20,agn,ge);
% {overlap compartment}
dx(21,agn,ge) = da(5,agn,ge)*x(5,agn,ge)+nu(agn,ge)*x(18,agn,ge) +chi(agn,ge)*x(19,agn,ge)+ii(agn,ge)*x(24,agn,ge)- lambda(21,agn,ge)*x(21,agn,ge)- epsilon(agn,ge)*x(21,agn,ge)- gamma(agn,ge)*x(21,agn,ge)-kappa(agn,ge)*x(21,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(21,agn,ge)+eta(agn-1)*x(21,agn-1,ge)-eta(agn)*x(21,agn,ge);%+Immigration(agn,ge)*(x(5,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(22,agn,ge) = da(6,agn,ge)*x(6,agn,ge)+etaa(agn,ge)*x(18,agn,ge) +psi(agn,ge)*x(20,agn,ge)+iota(agn,ge)*x(24,agn,ge)- lambda(22,agn,ge)*x(22,agn,ge)-Leps(agn,ge)*x(22,agn,ge)-theta(agn,ge)*x(22,agn,ge)-sampi(agn,ge)*x(22,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(22,agn,ge)+eta(agn-1)*x(22,agn-1,ge)-eta(agn)*x(22,agn,ge);%+Immigration(agn,ge)*(x(6,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(23,agn,ge) = da(7,agn,ge)*x(7,agn,ge)+omega(agn,ge)*x(19,agn,ge) +xi(agn,ge)*x(20,agn,ge)+upsilon(agn,ge)*x(24,agn,ge)- lambda(23,agn,ge)*x(23,agn,ge)- pii(agn,ge)*x(23,agn,ge)- rho(agn,ge)*x(23,agn,ge)-Omg(agn,ge)*x(23,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(23,agn,ge)+eta(agn-1)*x(23,agn-1,ge)-eta(agn)*x(23,agn,ge);%+Immigration(agn,ge)*(x(7,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(24,agn,ge) = da(8,agn,ge)*x(8,agn,ge)+kappa(agn,ge)*x(21,agn,ge)+sampi(agn,ge)*x(22,agn,ge)+Omg(agn,ge)*x(23,agn,ge)- lambda(24,agn,ge)*x(24,agn,ge)- iota(agn,ge)*x(24,agn,ge)- ii(agn,ge)*x(24,agn,ge)-upsilon(agn,ge)*x(24,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(24,agn,ge)+eta(agn-1)*x(24,agn-1,ge)-eta(agn)*x(24,agn,ge);%+Immigration(agn,ge)*(x(8,agn,ge)/squeeze(sum(x(:,agn,ge),1)));

% {Infected DM}
dx(25,agn,ge) =lambda(17,agn,ge)*x(17,agn,ge)+r(agn,ge)*x(26,agn,ge) + q(agn,ge)*x(27,agn,ge)+p(agn,ge)*x(28,agn,ge)- m(agn,ge)*x(25,agn,ge) - n(agn,ge)*x(25,agn,ge)-oo(agn,ge)*x(25,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(25,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(25,agn,ge)+eta(agn-1)*x(25,agn-1,ge)-eta(agn)*x(25,agn,ge);
dx(26,agn,ge) = lambda(18,agn,ge)*x(18,agn,ge) + a(agn,ge)*x(29,agn,ge) + b(agn,ge)*x(30,agn,ge)+m(agn,ge)*x(25,agn,ge)- e(agn,ge)*x(26,agn,ge) - ff(agn,ge)*x(26,agn,ge)-r(agn,ge)*x(26,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(26,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(26,agn,ge)+eta(agn-1)*x(26,agn-1,ge)-eta(agn)*x(26,agn,ge);
dx(27,agn,ge)= lambda(19,agn,ge)*x(19,agn,ge) + n(agn,ge)*x(25,agn,ge)+iii(agn,ge)*x(29,agn,ge)+cc(agn,ge)*x(31,agn,ge) - jj(agn,ge)*x(27,agn,ge) - gg(agn,ge)*x(27,agn,ge) -q(agn,ge)*x(27,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(27,agn,ge) -(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(27,agn,ge)+eta(agn-1)*x(27,agn-1,ge)-eta(agn)*x(27,agn,ge);
dx(28,agn,ge)=lambda(20,agn,ge)*x(20,agn,ge) + oo(agn,ge)*x(25,agn,ge)+k(agn,ge)*x(30,agn,ge)+ d(agn,ge)*x(31,agn,ge) - p(agn,ge)*x(28,agn,ge) - l(agn,ge)*x(28,agn,ge) - hh(agn,ge)*x(28,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(28,agn,ge) -(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(28,agn,ge)+eta(agn-1)*x(28,agn-1,ge)-eta(agn)*x(28,agn,ge);

dx(29,agn,ge)=lambda(21,agn,ge)*x(21,agn,ge) + ss(agn,ge)*x(32,agn,ge)+ jj(agn,ge)*x(27,agn,ge) + e(agn,ge)*x(26,agn,ge) - xx(agn,ge)*x(29,agn,ge) - a(agn,ge)*x(29,agn,ge)-iii(agn,ge)*x(29,agn,ge)... 
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(29,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(29,agn,ge)+eta(agn-1)*x(29,agn-1,ge)-eta(agn)*x(29,agn,ge);
dx(30,agn,ge)=lambda(22,agn,ge)*x(22,agn,ge) + tt(agn,ge)*x(32,agn,ge)+l(agn,ge)*x(28,agn,ge) + ff(agn,ge)*x(26,agn,ge)- w(agn,ge)*x(30,agn,ge)-b(agn,ge)*x(30,agn,ge)-k(agn,ge)*x(30,agn,ge)... 
- Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(30,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(30,agn,ge)+eta(agn-1)*x(30,agn-1,ge)-eta(agn)*x(30,agn,ge);
dx(31,agn,ge)=lambda(23,agn,ge)*x(23,agn,ge) + u(agn,ge)*x(32,agn,ge)+hh(agn,ge)*x(28,agn,ge) + gg(agn,ge)*x(27,agn,ge)- d(agn,ge)*x(31,agn,ge)-cc(agn,ge)*x(31,agn,ge)-v(agn,ge)*x(31,agn,ge)... 
- Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(31,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(31,agn,ge)+eta(agn-1)*x(31,agn-1,ge)-eta(agn)*x(31,agn,ge);

dx(32,agn,ge)=lambda(24,agn,ge)*x(24,agn,ge) +xx(agn,ge)*x(29,agn,ge) + w(agn,ge)*x(30,agn,ge)+v(agn,ge)*x(31,agn,ge) - ss(agn,ge)*x(32,agn,ge) -  tt(agn,ge)*x(32,agn,ge)- u(agn,ge)*x(32,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(32,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(32,agn,ge)+eta(agn-1)*x(32,agn-1,ge)-eta(agn)*x(32,agn,ge);
for gr=1:32
dmmort(1,agn,ge,gr)=Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(gr,agn,ge);
end
% dmmort(1,agn,ge)=(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge)+...
%                  (RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(25,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(26,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(27,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(28,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(29,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(30,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(31,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(32,agn,ge);
% dmmortYLL(1,agn,ge)=(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge);

 end
     end
else
    
dx(1,1,1) =N(1,1).*BirthFuncGaussian(t,a1f,b1f,c1f)-Mortality(t,1,b2f(1,1),c2f(1,1),a3f(1,1),b3f(1,1),c3f(1,1),s(1,1))*x(1,1,1)-eta(1)*x(1,1,1); %*brate*N*mu(agn,ge)*N0
dx(1,1,2) =N(2,1).*BirthFuncGaussian(t,a1m,b1m,c1m)-Mortality(t,1,b2f(1,2),c2f(1,2),a3f(1,2),b3f(1,2),c3f(1,2),s(1,2))*x(1,1,2)-eta(1)*x(1,1,2); %*brate*N*mu(agn,ge)*N0   
% dmmort(1,1)=0;
% dmmort(1,2)=0;
    for agn=2:ag
 for ge=1:gen
dx(1,agn,ge) = -da(1,agn,ge)*x(1,agn,ge)+sigma(agn,ge)*x(2,agn,ge)+delta(agn,ge)*x(3,agn,ge)+phi(agn,ge)*x(4,agn,ge)-lambda(1,agn,ge)*x(1,agn,ge)- beta(agn,ge)*x(1,agn,ge)-alpha(agn,ge)*x(1,agn,ge)-shei(agn,ge)*x(1,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(1,agn,ge)+eta(agn-1)*x(1,agn-1,ge)-eta(agn)*x(1,agn,ge);%+Immigration(agn,ge)*(x(1,agn,ge)/squeeze(sum(x(:,agn,ge),1))); 
% {Obese}
dx(2,agn,ge) = -da(2,agn,ge)*x(2,agn,ge)+alpha(agn,ge)*x(1,agn,ge)+ epsilon(agn,ge)*x(5,agn,ge)+theta(agn,ge)*x(6,agn,ge)-lambda(2,agn,ge)*x(2,agn,ge)-nu(agn,ge)*x(2,agn,ge)- etaa(agn,ge)*x(2,agn,ge)-sigma(agn,ge)*x(2,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(2,agn,ge)+eta(agn-1)*x(2,agn-1,ge)-eta(agn)*x(2,agn,ge);%+Immigration(agn,ge)*(x(2,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
% {Smoker}
dx(3,agn,ge) = -da(3,agn,ge)*x(3,agn,ge)+beta(agn,ge)*x(1,agn,ge)+gamma(agn,ge)*x(5,agn,ge)+pii(agn,ge)*x(7,agn,ge)-lambda(3,agn,ge)*x(3,agn,ge)- chi(agn,ge)*x(3,agn,ge) - delta(agn,ge)*x(3,agn,ge)-omega(agn,ge)*x(3,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(3,agn,ge)+eta(agn-1)*x(3,agn-1,ge)-eta(agn)*x(3,agn,ge);%+Immigration(agn,ge)*(x(3,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
% {Physical Active};
dx(4,agn,ge) =  -da(4,agn,ge)*x(4,agn,ge)+shei(agn,ge)*x(1,agn,ge)+rho(agn,ge)*x(7,agn,ge)+Leps(agn,ge)*x(6,agn,ge)-lambda(4,agn,ge)*x(4,agn,ge)  - phi(agn,ge)*x(4,agn,ge) - psi(agn,ge)*x(4,agn,ge)-xi(agn,ge)*x(4,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(4,agn,ge)+eta(agn-1)*x(4,agn-1,ge)-eta(agn)*x(4,agn,ge);%+Immigration(agn,ge)*(x(4,agn,ge)/squeeze(sum(x(:,agn,ge),1)));

% {overlap compartment}
dx(5,agn,ge) = - da(5,agn,ge)*x(5,agn,ge)+nu(agn,ge)*x(2,agn,ge) +chi(agn,ge)*x(3,agn,ge)+ii(agn,ge)*x(8,agn,ge)- lambda(5,agn,ge)*x(5,agn,ge)- epsilon(agn,ge)*x(5,agn,ge)- gamma(agn,ge)*x(5,agn,ge)-kappa(agn,ge)*x(5,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(5,agn,ge)+eta(agn-1)*x(5,agn-1,ge)-eta(agn)*x(5,agn,ge);%+Immigration(agn,ge)*(x(5,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(6,agn,ge) = - da(6,agn,ge)*x(6,agn,ge)+etaa(agn,ge)*x(2,agn,ge) +psi(agn,ge)*x(4,agn,ge)+iota(agn,ge)*x(8,agn,ge)- lambda(6,agn,ge)*x(6,agn,ge)-Leps(agn,ge)*x(6,agn,ge)-theta(agn,ge)*x(6,agn,ge)-sampi(agn,ge)*x(6,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(6,agn,ge)+eta(agn-1)*x(6,agn-1,ge)-eta(agn)*x(6,agn,ge);%+Immigration(agn,ge)*(x(6,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(7,agn,ge) =- da(7,agn,ge)*x(7,agn,ge)+ omega(agn,ge)*x(3,agn,ge) +xi(agn,ge)*x(4,agn,ge)+upsilon(agn,ge)*x(8,agn,ge)- lambda(7,agn,ge)*x(7,agn,ge)- pii(agn,ge)*x(7,agn,ge)- rho(agn,ge)*x(7,agn,ge)-Omg(agn,ge)*x(7,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(7,agn,ge)+eta(agn-1)*x(7,agn-1,ge)-eta(agn)*x(7,agn,ge);%+Immigration(agn,ge)*(x(7,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(8,agn,ge) =- da(8,agn,ge)*x(8,agn,ge)+ kappa(agn,ge)*x(5,agn,ge)+sampi(agn,ge)*x(6,agn,ge)+Omg(agn,ge)*x(7,agn,ge)- lambda(8,agn,ge)*x(8,agn,ge)- iota(agn,ge)*x(8,agn,ge)- ii(agn,ge)*x(8,agn,ge)-upsilon(agn,ge)*x(8,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(8,agn,ge)+eta(agn-1)*x(8,agn-1,ge)-eta(agn)*x(8,agn,ge);%+Immigration(agn,ge)*(x(8,agn,ge)/squeeze(sum(x(:,agn,ge),1)));

% {Infected DM}
dx(9,agn,ge) = lambda(1,agn,ge)*x(1,agn,ge) + r(agn,ge)*x(10,agn,ge) + q(agn,ge)*x(11,agn,ge)+p(agn,ge)*x(12,agn,ge)- m(agn,ge)*x(9,agn,ge) - n(agn,ge)*x(9,agn,ge)-oo(agn,ge)*x(9,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)+eta(agn-1)*x(9,agn-1,ge)...
    -eta(agn)*x(9,agn,ge);%+Immigration(agn,ge)*(x(9,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(10,agn,ge) = lambda(2,agn,ge)*x(2,agn,ge) + a(agn,ge)*x(13,agn,ge) + b(agn,ge)*x(14,agn,ge)+m(agn,ge)*x(9,agn,ge)- e(agn,ge)*x(10,agn,ge) - ff(agn,ge)*x(10,agn,ge)-r(agn,ge)*x(10,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)+eta(agn-1)*x(10,agn-1,ge)...
    -eta(agn)*x(10,agn,ge);%+Immigration(agn,ge)*(x(10,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(11,agn,ge)= lambda(3,agn,ge)*x(3,agn,ge) + n(agn,ge)*x(9,agn,ge)+iii(agn,ge)*x(13,agn,ge)+cc(agn,ge)*x(15,agn,ge) - jj(agn,ge)*x(11,agn,ge) - gg(agn,ge)*x(11,agn,ge) -q(agn,ge)*x(11,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge) -(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge)+eta(agn-1)*x(11,agn-1,ge)...
    -eta(agn)*x(11,agn,ge);%+Immigration(agn,ge)*(x(11,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(12,agn,ge)= lambda(4,agn,ge)*x(4,agn,ge) + oo(agn,ge)*x(9,agn,ge)+k(agn,ge)*x(14,agn,ge)+ d(agn,ge)*x(15,agn,ge) - p(agn,ge)*x(12,agn,ge) - l(agn,ge)*x(12,agn,ge) - hh(agn,ge)*x(12,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge) -(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge)+eta(agn-1)*x(12,agn-1,ge)...
    -eta(agn)*x(12,agn,ge);%+Immigration(agn,ge)*(x(12,agn,ge)/squeeze(sum(x(:,agn,ge),1)));

dx(13,agn,ge)= lambda(5,agn,ge)*x(5,agn,ge) + ss(agn,ge)*x(16,agn,ge)+ jj(agn,ge)*x(11,agn,ge) + e(agn,ge)*x(10,agn,ge) - xx(agn,ge)*x(13,agn,ge) - a(agn,ge)*x(13,agn,ge)-iii(agn,ge)*x(13,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)+eta(agn-1)*x(13,agn-1,ge)...
    -eta(agn)*x(13,agn,ge);%+Immigration(agn,ge)*(x(13,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(14,agn,ge)= lambda(6,agn,ge)*x(6,agn,ge) + tt(agn,ge)*x(16,agn,ge)+l(agn,ge)*x(12,agn,ge) + ff(agn,ge)*x(10,agn,ge)- w(agn,ge)*x(14,agn,ge)-b(agn,ge)*x(14,agn,ge)-k(agn,ge)*x(14,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)+eta(agn-1)*x(14,agn-1,ge)...
    -eta(agn)*x(14,agn,ge);%+Immigration(agn,ge)*(x(14,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(15,agn,ge)= lambda(7,agn,ge)*x(7,agn,ge) + u(agn,ge)*x(16,agn,ge)+hh(agn,ge)*x(12,agn,ge) + gg(agn,ge)*x(11,agn,ge)- d(agn,ge)*x(15,agn,ge)-cc(agn,ge)*x(15,agn,ge)-v(agn,ge)*x(15,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)+eta(agn-1)*x(15,agn-1,ge)...
    -eta(agn)*x(15,agn,ge);%+Immigration(agn,ge)*(x(15,agn,ge)/squeeze(sum(x(:,agn,ge),1)));

dx(16,agn,ge)= lambda(8,agn,ge)*x(8,agn,ge) +xx(agn,ge)*x(13,agn,ge) + w(agn,ge)*x(14,agn,ge)+v(agn,ge)*x(15,agn,ge) - ss(agn,ge)*x(16,agn,ge) -  tt(agn,ge)*x(16,agn,ge)- u(agn,ge)*x(16,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge)+eta(agn-1)*x(16,agn-1,ge)...
    -eta(agn)*x(16,agn,ge);%+Immigration(agn,ge)*(x(16,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
%% invtervention
dx(17,agn,ge) = da(1,agn,ge)*x(1,agn,ge)+sigma(agn,ge)*x(18,agn,ge)+delta(agn,ge)*x(19,agn,ge)+phi(agn,ge)*x(20,agn,ge)-lambda(17,agn,ge)*x(17,agn,ge)- beta(agn,ge)*x(17,agn,ge)-alpha(agn,ge)*x(17,agn,ge)-shei(agn,ge)*x(17,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(17,agn,ge)+eta(agn-1)*x(17,agn-1,ge)-eta(agn)*x(17,agn,ge);%+Immigration(agn,ge)*(x(1,agn,ge)/squeeze(sum(x(:,agn,ge),1))); 
% {Obese}
dx(18,agn,ge) = da(2,agn,ge)*x(2,agn,ge)+alpha(agn,ge)*x(17,agn,ge)+ epsilon(agn,ge)*x(21,agn,ge)+theta(agn,ge)*x(22,agn,ge)-lambda(18,agn,ge)*x(18,agn,ge)-nu(agn,ge)*x(18,agn,ge)-etaa(agn,ge)*x(18,agn,ge)-sigma(agn,ge)*x(18,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(18,agn,ge)+eta(agn-1)*x(18,agn-1,ge)-eta(agn)*x(18,agn,ge);%+Immigration(agn,ge)*(x(agn,ge)/squeeze(sum(x(:,agn,ge),1)));
% {Smoker}
dx(19,agn,ge) = da(3,agn,ge)*x(3,agn,ge)+beta(agn,ge)*x(17,agn,ge)+gamma(agn,ge)*x(21,agn,ge)+pii(agn,ge)*x(23,agn,ge)-lambda(19,agn,ge)*x(19,agn,ge)- chi(agn,ge)*x(19,agn,ge)-delta(agn,ge)*x(19,agn,ge)-omega(agn,ge)*x(19,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(19,agn,ge)+eta(agn-1)*x(19,agn-1,ge)-eta(agn)*x(19,agn,ge);%+Immigration(agn,ge)*(x(3,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
% {Physical Active};
dx(20,agn,ge) = da(4,agn,ge)*x(4,agn,ge)+shei(agn,ge)*x(17,agn,ge)+rho(agn,ge)*x(23,agn,ge)+Leps(agn,ge)*x(22,agn,ge)-lambda(20,agn,ge)*x(20,agn,ge)- phi(agn,ge)*x(20,agn,ge)-psi(agn,ge)*x(20,agn,ge)-xi(agn,ge)*x(20,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(20,agn,ge)+eta(agn-1)*x(20,agn-1,ge)-eta(agn)*x(20,agn,ge);%+Immigration(agn,ge)*(x(4,agn,ge)/squeeze(sum(x(:,agn,ge),1)));

% {overlap compartment}
dx(21,agn,ge) = da(5,agn,ge)*x(5,agn,ge)+nu(agn,ge)*x(18,agn,ge) +chi(agn,ge)*x(19,agn,ge)+ii(agn,ge)*x(24,agn,ge)- lambda(21,agn,ge)*x(21,agn,ge)- epsilon(agn,ge)*x(21,agn,ge)- gamma(agn,ge)*x(21,agn,ge)-kappa(agn,ge)*x(21,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(21,agn,ge)+eta(agn-1)*x(21,agn-1,ge)-eta(agn)*x(21,agn,ge);%+Immigration(agn,ge)*(x(5,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(22,agn,ge) = da(6,agn,ge)*x(6,agn,ge)+etaa(agn,ge)*x(18,agn,ge) +psi(agn,ge)*x(20,agn,ge)+iota(agn,ge)*x(24,agn,ge)- lambda(22,agn,ge)*x(22,agn,ge)-Leps(agn,ge)*x(22,agn,ge)-theta(agn,ge)*x(22,agn,ge)-sampi(agn,ge)*x(22,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(22,agn,ge)+eta(agn-1)*x(22,agn-1,ge)-eta(agn)*x(22,agn,ge);%+Immigration(agn,ge)*(x(6,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(23,agn,ge) = da(7,agn,ge)*x(7,agn,ge)+omega(agn,ge)*x(19,agn,ge) +xi(agn,ge)*x(20,agn,ge)+upsilon(agn,ge)*x(24,agn,ge)- lambda(23,agn,ge)*x(23,agn,ge)- pii(agn,ge)*x(23,agn,ge)- rho(agn,ge)*x(23,agn,ge)-Omg(agn,ge)*x(23,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(23,agn,ge)+eta(agn-1)*x(23,agn-1,ge)-eta(agn)*x(23,agn,ge);%+Immigration(agn,ge)*(x(7,agn,ge)/squeeze(sum(x(:,agn,ge),1)));
dx(24,agn,ge) = da(8,agn,ge)*x(8,agn,ge)+kappa(agn,ge)*x(21,agn,ge)+sampi(agn,ge)*x(22,agn,ge)+Omg(agn,ge)*x(23,agn,ge)- lambda(24,agn,ge)*x(24,agn,ge)- iota(agn,ge)*x(24,agn,ge)- ii(agn,ge)*x(24,agn,ge)-upsilon(agn,ge)*x(24,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(24,agn,ge)+eta(agn-1)*x(24,agn-1,ge)-eta(agn)*x(24,agn,ge);%+Immigration(agn,ge)*(x(8,agn,ge)/squeeze(sum(x(:,agn,ge),1)));

% {Infected DM}

dx(25,agn,ge) =lambda(17,agn,ge)*x(17,agn,ge)+r(agn,ge)*x(26,agn,ge) + q(agn,ge)*x(27,agn,ge)+p(agn,ge)*x(28,agn,ge)- m(agn,ge)*x(25,agn,ge) - n(agn,ge)*x(25,agn,ge)-oo(agn,ge)*x(25,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(25,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(25,agn,ge)+eta(agn-1)*x(25,agn-1,ge)-eta(agn)*x(25,agn,ge);
dx(26,agn,ge) = lambda(18,agn,ge)*x(18,agn,ge) + a(agn,ge)*x(29,agn,ge) + b(agn,ge)*x(30,agn,ge)+m(agn,ge)*x(25,agn,ge)- e(agn,ge)*x(26,agn,ge) - ff(agn,ge)*x(26,agn,ge)-r(agn,ge)*x(26,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(26,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(26,agn,ge)+eta(agn-1)*x(26,agn-1,ge)-eta(agn)*x(26,agn,ge);
dx(27,agn,ge)= lambda(19,agn,ge)*x(19,agn,ge) + n(agn,ge)*x(25,agn,ge)+iii(agn,ge)*x(29,agn,ge)+cc(agn,ge)*x(31,agn,ge) - jj(agn,ge)*x(27,agn,ge) - gg(agn,ge)*x(27,agn,ge) -q(agn,ge)*x(27,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(27,agn,ge) -(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(27,agn,ge)+eta(agn-1)*x(27,agn-1,ge)-eta(agn)*x(27,agn,ge);
dx(28,agn,ge)=lambda(20,agn,ge)*x(20,agn,ge) + oo(agn,ge)*x(25,agn,ge)+k(agn,ge)*x(30,agn,ge)+ d(agn,ge)*x(31,agn,ge) - p(agn,ge)*x(28,agn,ge) - l(agn,ge)*x(28,agn,ge) - hh(agn,ge)*x(28,agn,ge)...
    -Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(28,agn,ge) -(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(28,agn,ge)+eta(agn-1)*x(28,agn-1,ge)-eta(agn)*x(28,agn,ge);

dx(29,agn,ge)=lambda(21,agn,ge)*x(21,agn,ge) + ss(agn,ge)*x(32,agn,ge)+ jj(agn,ge)*x(27,agn,ge) + e(agn,ge)*x(26,agn,ge) - xx(agn,ge)*x(29,agn,ge) - a(agn,ge)*x(29,agn,ge)-iii(agn,ge)*x(29,agn,ge)... 
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(29,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(29,agn,ge)+eta(agn-1)*x(29,agn-1,ge)-eta(agn)*x(29,agn,ge);
dx(30,agn,ge)=lambda(22,agn,ge)*x(22,agn,ge) + tt(agn,ge)*x(32,agn,ge)+l(agn,ge)*x(28,agn,ge) + ff(agn,ge)*x(26,agn,ge)- w(agn,ge)*x(30,agn,ge)-b(agn,ge)*x(30,agn,ge)-k(agn,ge)*x(30,agn,ge)... 
- Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(30,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(30,agn,ge)+eta(agn-1)*x(30,agn-1,ge)-eta(agn)*x(30,agn,ge);
dx(31,agn,ge)=lambda(23,agn,ge)*x(23,agn,ge) + u(agn,ge)*x(32,agn,ge)+hh(agn,ge)*x(28,agn,ge) + gg(agn,ge)*x(27,agn,ge)- d(agn,ge)*x(31,agn,ge)-cc(agn,ge)*x(31,agn,ge)-v(agn,ge)*x(31,agn,ge)... 
- Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(31,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(31,agn,ge)+eta(agn-1)*x(31,agn-1,ge)-eta(agn)*x(31,agn,ge);

dx(32,agn,ge)=lambda(24,agn,ge)*x(24,agn,ge) +xx(agn,ge)*x(29,agn,ge) + w(agn,ge)*x(30,agn,ge)+v(agn,ge)*x(31,agn,ge) - ss(agn,ge)*x(32,agn,ge) -  tt(agn,ge)*x(32,agn,ge)- u(agn,ge)*x(32,agn,ge)...
    - Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(32,agn,ge)-(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(32,agn,ge)+eta(agn-1)*x(32,agn-1,ge)-eta(agn)*x(32,agn,ge);

% dmmort(1,agn,ge)=(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge)+...
                 (RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(25,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(26,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(27,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(28,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(29,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(30,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(31,agn,ge)+(RRm(agn,ge))*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(32,agn,ge);
% dmmortYLL(1,agn,ge)=(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge);
for gr=1:32
dmmort(1,agn,ge,gr)=Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(gr,agn,ge);
end

 end
    end
end
% for agn=1:ag
%     for ge=1:gen
%         dmmort(1,agn,ge)=(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge)+...
%                  (RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(25,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(26,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(27,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(28,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(29,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(30,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(31,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(32,agn,ge);
% dmmortYLL(1,agn,ge)=(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(9,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(10,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(11,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(12,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(13,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(14,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(15,agn,ge)+(RRm(agn,ge)-1)*Mortality(t,agn,b2f(1,ge),c2f(1,ge),a3f(1,ge),b3f(1,ge),c3f(1,ge),s(1,ge))*x(16,agn,ge);
%     end
% end
dxL=reshape(dx,numel(dx),1);
end