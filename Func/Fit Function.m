% Function to fit a ACD(q,p) model to data

function [specOut]=ACD_Fit(x,dist,q,p,method,indep)

if size(x,2)>1
    error('The input x should be a vector')
end

if strcmp(dist,'exp')==0&&strcmp(dist,'weibull')==0
    error('The input dist should be either ''exp'' or ''weibull''');
end

if nargin()==4
    method=1;
    indep=[];
end

if nargin()==5
    indep=[];
end

if ~isempty(indep)
    if size(indep,1)~=size(x,1)
        error('The number of rows in x should match the number of rows in indep');
    end
end

if ~any(method==[1,2,3,4])
    error(' The input method should be 1,2,3 or 4.');
end

if q<1||p<1
    error(' The input q and p should be integers higher or equal than one');
end

% Some precalculation for param0

nIndep=size(indep,2);

for i=0:q-1
    OLS_indep(:,i+1)=x(1+i:end-q+i);
end

param_OLS=regress(x(q+1:end,1),[ones(length(x)-q,1) , OLS_indep,indep(max(p,q)+1:end,:)]); % simple OLS for alpha0 and beta0

switch dist

    case 'exp'

        if ~isempty(indep)
            param0=[.05 param_OLS(2:1+q)' repmat((1-sum(param_OLS(2:1+q)))/p,1,p) param_OLS(1+q+1:end)'];
        else
            param0=[.05 param_OLS(2:end)' repmat((1-sum(param_OLS(2:1+q)))/p,1,p)];
        end

    case 'weibull'

        if ~isempty(indep)
            param0=[.05 param_OLS(2:1+q)' repmat((1-sum(param_OLS(2:1+q)))/p,1,p) param_OLS(1+q+1:end)' .8];
        else
            param0=[.05 param_OLS(2:end)' repmat((1-sum(param_OLS(2:1+q)))/p,1,p) .8];
        end


end

options=optimset('fminsearch');
options=optimset(options,'display','off');

warning('off');

global global_p; 
global global_q;

global_p=p;
global_q=q;

[param]=fminsearch(@(param)ACD_Lik(param,x,q,p,dist,1,indep),param0,options);

[V]=getvarMatrix_ACD(method,param,x,q,p,dist,indep);

param_std=sqrt(diag(V));

Coeff.w_std=param_std(1);
Coeff.q_std=param_std(2:q+1);
Coeff.p_std=param_std(q+2:q+1+p);

Coeff.indep_c_std=param_std(q+1+p+1:q+1+p+nIndep);

if strmatch(dist,'weibull')
    Coeff.y_std=param_std(end);
end

[sumLik,specOut]=ACD_Lik(param,x,q,p,dist,0,indep); % filtering it again in order to recover the conditional durations

fprintf(1,'\n\n******* Optimization Finished *******\n\n');
fprintf(1,['Maximum Log Likelihood: ' num2str(-sumLik) '\n\n']);

fprintf(1,['Parameters for ACD(' num2str(q) ',' num2str(p) ') Model (std error in parenthesis):'])
fprintf(1,['\n  Const (Coeff.w) = ' num2str(specOut.w,'%2.4f') ' (' num2str(Coeff.w_std,'%2.4f') ')']);

for i=1:q
    fprintf(1,['\n  Alpha ' num2str(i) ' (Coeff.q(' num2str(i) ')) = ' num2str(specOut.q(i),'%2.4f') ' (' num2str(Coeff.q_std(i),'%2.4f') ')']);
end
for i=1:p
    fprintf(1,['\n  Beta  ' num2str(i) ' (Coeff.p(' num2str(i) ')) = ' num2str(specOut.p(i),'%2.4f') ' (' num2str(Coeff.p_std(i),'%2.4f') ')']);
end

if ~isempty(indep)
    for i=1:nIndep
        fprintf(1,['\n  Indep Col  ' num2str(i) ' (Coeff.beta_c(' num2str(i) ')) = ' num2str(specOut.indep_c(i),'%2.4f') ' (' num2str(Coeff.indep_c_std(i),'%2.4f') ')']);
    end
end


switch dist
    case 'weibull'
        fprintf(1,['\n  Weibull param (Coeff.y) = ' num2str(specOut.y,'%2.4f') '(' num2str(Coeff.y_std,'%2.4f') ')']);
        fprintf(1,'\n\n');
    case 'exp'
        fprintf(1,'\n\n');
end
