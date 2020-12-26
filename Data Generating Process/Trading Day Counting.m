% Import HF Data.xls as matrix into Matlab %
% The code is to find out the number of trading days in the dataset %

i = 1;
for n = 1:size(HFData,1)-1
    
    if HFData(n+1,1) ~= HFData(n,1)
        i = i+1;
    end
    
end

fprintf('%d',i);
