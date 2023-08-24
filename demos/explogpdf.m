function v = explogpdf(x,beta,p)
arguments
    x (:,1) double
    beta (1,1) double
    p (1,1) double {mustBeGreaterThan(p,1)} = 1+eps
end 

v = (-1./log(p)).* beta.*(1-p).*exp(-beta.*x)./(1-(1-p).*exp(-beta.*x));
