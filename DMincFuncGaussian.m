function [lambdaH] =DMincFuncGaussian(ag,gen,a,b,c)
lambdass=zeros(ag,gen);
for age=1:ag 
    for ge=1:gen
lambdass(age,ge)=a(1,ge)*exp(-((age-b(1,ge))/c(1,ge))^2);
    end
end%
lambdaH=lambdass';
end
    


