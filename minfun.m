function Cost = minfun(par_t,Final_param1,data1, BFf, BMm,Btot)
% global BFf BMm Btot
% Transform from whole real range to sensible ranges for parameters:
% scale = [1e-4 0.2 5 10 100 30 1];
fitpara_t = par_t;%.*scale';
% params_dummy=[];
% params_dummy.rskpar = zeros(1,length(fitpara_t));
% for i = 1:length(fitpara_t)
% params_dummy.rskpar(i) = erfb(ranges(i,1),ranges(i,2),fitpara_t(i));
% end
% w=[ones(1,15),repmat(2,1,12),ones(1,76)];
[propQatar101F,propQatar102F,propQatar103F,propQatar104F,propQatar105F,propQatar106F,...  
    propQatar121F,propQatar125F,...
    propQatar131F,propQatar132F,propQatar133F,propQatar134F,propQatar135F,...
    propQatar151F,propQatar152F,propQatar155F,...
    propQatar190F,propQatar191F,propQatar192F,propQatar193F,propQatar194F,propQatar195F,propQatar196F,propQatar197F,propQatar198F,...
    propQatar101M,propQatar102M,propQatar103M,propQatar104M,propQatar105M,propQatar106M,...  
    propQatar121M,propQatar125M,...
    propQatar131M,propQatar132M,propQatar133M,propQatar134M,propQatar135M,...
    propQatar151M,propQatar152M,propQatar155M,...
    propQatar190M,propQatar191M,propQatar192M,propQatar193M,propQatar194M,propQatar195M,propQatar196M,propQatar197M,propQatar198M,...
    SumTOT,S1sumTotF,S1sumTotM]=DM_model_3RFsDEMO(fitpara_t,Final_param1,Btot);



% %%Least square
% cc=0;
% differ=0;
cc1F=0;
% differ1=0;
cc2F=0;
cc1M=0;
cc2M=0;


differ10F=((propQatar101F-BFf(1,1))/BFf(1,1))^2+((propQatar102F-BFf(1,2))/BFf(1,2))^2+((propQatar103F-BFf(1,3))/BFf(1,3))^2+...
    ((propQatar104F-BFf(1,4))/BFf(1,4))^2+((propQatar105F-BFf(1,5))/BFf(1,5))^2+((propQatar106F-BFf(1,6))/BFf(1,6))^2;
differ13F=((propQatar131F-BFf(3,1))/BFf(3,1))^2+((propQatar132F-BFf(3,2))/BFf(3,2))^2+((propQatar133F-BFf(3,3))/BFf(3,3))^2+...
    ((propQatar134F-BFf(3,4))/BFf(3,4))^2+((propQatar135F-BFf(3,5))/BFf(3,5))^2;
differ19F=((propQatar190F-BFf(5,2))/BFf(5,2))^2+((propQatar191F-BFf(5,3))/BFf(5,3))^2+((propQatar192F-BFf(5,4))/BFf(5,4))^2+...
    ((propQatar193F-BFf(5,5))/BFf(5,5))^2+((propQatar194F-BFf(5,6))/BFf(5,6))^2+...
    ((propQatar195F-BFf(5,7))/BFf(5,7))^2+((propQatar196F-BFf(5,8))/BFf(5,8))^2+...
    ((propQatar197F-BFf(5,9))/BFf(5,9))^2+((propQatar198F-BFf(5,10))/BFf(5,10))^2;

differ10M=((propQatar101M-BMm(1,1))/BMm(1,1))^2+((propQatar102M-BMm(1,2))/BMm(1,2))^2+((propQatar103M-BMm(1,3))/BMm(1,3))^2+...
    ((propQatar104M-BMm(1,4))/BMm(1,4))^2+((propQatar105M-BMm(1,5))/BMm(1,5))^2+((propQatar106M-BMm(1,6))/BMm(1,6))^2;
differ13M=((propQatar131M-BMm(3,1))/BMm(3,1))^2+((propQatar132M-BMm(3,2))/BMm(3,2))^2+((propQatar133M-BMm(3,3))/BMm(3,3))^2+...
    ((propQatar134M-BMm(3,4))/BMm(3,4))^2+((propQatar135M-BMm(3,5))/BMm(3,5))^2;
differ19M=((propQatar190M-BMm(5,2))/BMm(5,2))^2+((propQatar191M-BMm(5,3))/BMm(5,3))^2+((propQatar192M-BMm(5,4))/BMm(5,4))^2+...
    ((propQatar193M-BMm(5,5))/BMm(5,5))^2+((propQatar194M-BMm(5,6))/BMm(5,6))^2+...
    ((propQatar195M-BMm(5,7))/BMm(5,7))^2+((propQatar196M-BMm(5,8))/BMm(5,8))^2+...
    ((propQatar197M-BMm(5,9))/BMm(5,9))^2+((propQatar198M-BMm(5,10))/BMm(5,10))^2;

ccf12=0;
t0=1904;
Nage=20;    
    for a=1:10
        differ12F=ccf12+((propQatar121F(a)-BFf(2,a))/(BFf(2,a)))^2;
        ccf12=differ12F+((propQatar125F-BFf(2,11))/BFf(2,11))^2;
    end
    ccf15=0;
    for a=1:12
        differ15F=ccf15+((propQatar152F(a)-BFf(4,a+1))/(BFf(4,a+1)))^2;
        ccf15=differ15F;
    end
    ccf15=differ15F+((propQatar155F-BFf(4,14))/BFf(4,14))^2+((propQatar151F-BFf(4,1))/BFf(4,1))^2;
ccm12=0;
        for a=1:10
        differ12M=ccm12+((propQatar121M(a)-BFf(2,a))/(BMm(2,a)))^2;
        ccm12=differ12M+((propQatar125M-BMm(2,11))/BMm(2,11))^2;
        end
    
    ccm15=0;
    for age=1:12
        differ15M=ccm15+((propQatar152M(age)-BMm(4,age+1))/(BMm(4,age+1)))^2;
        ccm15=differ15M;
    end
    ccm15=differ15M+((propQatar155M-BMm(4,14))/BMm(4,14))^2+((propQatar151M-BMm(4,1))/BMm(4,1))^2;
% end
cc1=0;
cc2=0;
cc3=0;
for t=2020:1:2049
        differ1=cc1+(((SumTOT(t-t0+1,1))-Btot((t-1950)+1,3))/Btot((t-1950)+1,3))^2;
        cc1=differ1;
        differ2=cc2+(((S1sumTotM(t-t0+1,1))-Btot((t-1950)+1,2))/Btot((t-1950)+1,2))^2;
        cc2=differ2;
        differ3=cc3+(((S1sumTotF(t-t0+1,1))-Btot((t-1950)+1,1))/Btot((t-1950)+1,1))^2;
        cc3=differ3;
end
Cost=2*(differ2+differ3)+cc1+ccm15+ccm12+ccf15+ccf12+differ10F+differ13F+2*differ19F+differ10M+differ13M+4*differ19M;

% Error_Vector=Model_Curve-data1;
% sse=sum((Error_Vector.^2));%w.*
% dev=sse;