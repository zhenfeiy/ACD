% Likelihood function for ACD_Fit.m

function [sumLik,specOut,loglik]=ACD_Lik(param,x,q,p,dist,dpl,indep)

if nargin==5
    dpl=1;
    indep=[];
end

nIndep=size(indep,2);

% Organizing Coefficients

Coeff.w=param(1);
Coeff.q=param(2:q+1);
Coeff.p=param(q+2:q+1+p);

if ~isempty(indep)
    Coeff.indep_c=param(q+p+2:q+p+1+nIndep);
end

nr=size(x,1);

switch dist
    case 'weibull'
        Coeff.y=param(end); % the extra parameter of weibull distribution
end

if isempty(indep)   % if no indep matrix is found, use vectorized filter function
    h1=Coeff.w+filter(Coeff.q,1,x);     % Vectorized Filter for conditional Durations
    h=filter([0 1],[1 -Coeff.p],h1);    % Vectorized Filter for conditional Durations
else
    h=zeros(nr,1);
    h(1:max(p,q))=0;
    for i=max(p,q)+1:nr
        h(i,1)=Coeff.w + flipdim(Coeff.q,2)*(x(i-q:i-1,1))+flipdim(Coeff.p,2)*h(i-p:i-1,1) + Coeff.indep_c*indep(i,:)';   % this is the psi equation
    end

end

switch dist
    case 'exp'
        loglik(:,1)=log(1./h(:,1).*exp(-x(:,1)./h(:,1))); % this is the log likelihood
    case 'weibull'
        % this is the log likelihood
        loglik(:,1)=log(Coeff.y./x(:,1).*((x(:,1).*gamma(1+1/Coeff.y))./h(:,1)).^Coeff.y.*exp(-(((x(:,1).*(gamma(1+1/Coeff.y))./h(:,1)).^Coeff.y))));
end

% control for inf and nan

loglik(1:max(p,q))=[];

infIdx=isinf(loglik);
loglik(infIdx)=-inf;

nanIdx=isnan(loglik);
loglik(nanIdx)=-inf;

% fmincon minimizes it, so I need the negative log likelihood

sumLik=-sum(loglik);

if isnan(sumLik)||isreal(sumLik)==0||isinf(sumLik)
    sumLik=inf;
end

% building specOut structure

specOut.h=h;
specOut.w=param(1);
specOut.q=param(2:q+1);
specOut.p=param(q+2:q+1+p);

if ~isempty(indep)
    specOut.indep_c=param(q+1+p+1:q+p+1+nIndep);
end

switch dist
    case 'weibull'
        specOut.y=param(end);
end

if dpl
    fprintf(1,['Log Likelihood ' dist ' ACD=%4.4f\n'],-sumLik);
end
