function prob_default = prob_default(I,m,p,rho,flag)
% flag = 1 -> nchoosek
% flag = 2 -> Stirling

M = 6; % integration bound

K = norminv(p);
p_y = @(y) normcdf((K-sqrt(rho)*y)/sqrt(1-rho));

if flag == 1
    integrand_HP = @(y) normpdf(y).*nchoosek(I,m).*(p_y(y).^m).*(1-p_y(y)).^(I-m); 
    prob_default = quadgk(integrand_HP,-M,M);
elseif flag == 2
    
    
end


end