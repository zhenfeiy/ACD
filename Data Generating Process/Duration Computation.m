load('HF_datetime');
load('HF_price');

i = 1;
Duration = [];
Dur_datetime = string([]);
temp = HF_datetime(1);

for j = 1:size(HF_datetime,1)
    if ~isnan(HF_price(j))
        value = HF_price(j);
    else
        HF_price(j) = value;
    end
end

for n = 1:size(HF_datetime,1)-1
    t1_test = datevec(datenum(HF_datetime(n)));
    t2_test = datevec(datenum(HF_datetime(n+1)));
    time_interval_test = etime(t2_test,t1_test)/60;
    
    if abs(time_interval_test) >= 90
        temp = HF_datetime(n+1);
        continue
    end
        
    if abs(HF_price(n+1)-HF_price(n)) >= 0.01
        t1 = datevec(datenum(temp));
        t2 = datevec(datenum(HF_datetime(n+1)));
        
        time_interval = etime(t2,t1)/60;
                
        Duration(i,1) = time_interval;
        Dur_datetime(i,1) = HF_datetime(n+1);
        i = i+1;
        temp = HF_datetime(n+1);
    end   
end
plot(Duration);
save('Dur_datetime.mat');
save('Duration.mat');

xlswrite('Dur_datetime.xlsx',Dur_datetime);
xlswrite('Duration.xlsx',Duration);
