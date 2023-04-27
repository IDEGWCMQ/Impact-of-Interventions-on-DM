function Cost = DMminfun(par_t,Result,data1, BFf, BMm,Btot)

fitpara_t = par_t;%.*scale';

[Model_Curve]=DM_model_3RFsDMQatar(fitpara_t,Result, BFf, BMm,Btot);

cc=0;
%% Least square: 
% Error_Vector=Model_Curve-data1;
% sse=sum((Error_Vector.^2));%w.*
% dev=sse;
for vt=1:length(data1)
  differ=cc+((Model_Curve(vt,1)-data1(vt,1))./data1(vt,1))^2;
  cc=differ;             
 end
   
Cost=differ;
save('fittingDMparam','fitpara_t')
% 
