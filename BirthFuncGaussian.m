%function [gammadi] =BirthFuncGaussian(t,Mort,a,b,c,MortBiss)
function [gammadi] =BirthFuncGaussian(t,a,b,c)
    gammadi=a*exp(-((t-b)/c)^2);
    %gammadi=a;
    


