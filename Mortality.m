function [gammadi] =Mortality(t,aa,b2,c2,a3,b3,c3,s) 
% if t<1985
        gammadi=(a3*exp(-((t-b3)/c3)^2))*(1/(s*(1+exp(-c2*(aa-b2))))); 
% else 
%             gammadi=s*(1/((1+exp(-c2*(aa-b2))))); 
%    gammadi=(a1*exp(-((t-b1)/c1)^2)); 
end
%gammadi=a1;