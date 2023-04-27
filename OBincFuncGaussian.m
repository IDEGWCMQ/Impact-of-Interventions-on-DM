%function [gammadi] =BirthFuncGaussian(t,Mort,a,b,c,MortBiss)
function [alpha] =OBincFuncGaussian(ag,gen,a,b,c)
for age=1:ag 
    for ge=1:gen
alpha(age,ge)=a(1,ge)*exp(-((age-b(1,ge))/c(1,ge))^2);
    end
end%gammadi=a;
    


