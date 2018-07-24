SSdat = [];
n = 100;
kon = [.1 .09 .08 .07 .06 .05 .04 .03 .02 .01];
koff = [.1 .2 .5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000];
for konloop =  [".1" ".09" ".08" ".07" ".06" ".05" ".04" ".03" ".02" ".01" ]
    for koffloop = [".1" ".2" ".5" "1" "2" "5" "10" "20" "50" "100" "200" "500" "1000" "2000" "5000" "10000" "20000" "50000" "100000"]
        file = "SmoldynDimerization/Data/Decay/rbind"+konloop+"/koff"+koffloop+"trim.txt";
        data = load(file);
        start = 1* floor(length(data)/3);
        monDat = data(start:end,1);
        SSdat = [SSdat,mean(monDat)/1000];
        clear data;
    end
end

SSdat = reshape(SSdat,[19,10]);

figure()
hold on
set(gca,'xscale','log')
xlabel('Koff values')
ylabel('Monomerization Fraction')
title('Steady State Monomerization Fractions in Smoldyn')
konEquiv = [];

for j = 1:10
    plot(koff,SSdat(:,j))
    konEquiv = [konEquiv,kon(j)];
    konEquiv = [konEquiv,interp1(SSdat(:,j),koff,.5,'pchip')];
end

legend(".1", ".09", ".08", ".07" ,".06", ".05", ".04", ".03", ".02", ".01" )