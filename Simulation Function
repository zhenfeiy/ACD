% Function to simulate a ACD(q,p) model to data

function [simulDur]=ACD_Simul(nr,Coeff,q,p,dist)

    if nr<(max(p,q)+1)
        error('The number of observations shouldbe be higher than max(q,p)')
    end
    
    
    if p<1||q<1
        error('The input q and p should be higher or equal than 1');
    end
    
    if strcmp(dist,'exp')==0&&strcmp(dist,'weibull')==0
        error('The input dist should be either ''exp'' or ''weibull''');
    end
    

    % Preallocation of large matrices and some precalculations
    
    simulDur=zeros(nr,1);
    durExp=zeros(nr,1);

    first_idx=max(q,p)+1;

    durExp(1,1)=0;
    simulDur(1:first_idx,1)=Coeff.w;
        
    for i=first_idx:nr

        % this is the psi equation
        
        durExp(i,1)=Coeff.w+flipdim(Coeff.q,2)*simulDur(i-q:i-1,1)+flipdim(Coeff.p,2)*durExp(i-p:i-1,1);
        
        switch dist
            
            case 'exp'
                simulDur(i,1)=durExp(i,1)*exprnd(1); % simulated duration at observation i
            case 'weibull'
                simulDur(i,1)=durExp(i,1)*wblrnd(1,Coeff.y);    % simulated duration at observation i

        end

    end

    plot(simulDur);
    title('Simulation of a Duration Series');
    xlabel('Observations');
    ylabel('Simulated Durations');
