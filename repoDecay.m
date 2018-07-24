function repoDecay(rbind,koff)
file = "SmoldynDimerization/Data/Decay/rbind"+rbind+"/koff"+koff+".txt";
data = load(file);
data = data/1000;
tdat = 1:length(data);
tdat = tdat * .00001;
figure()
plot(tdat,data(:,1))
xlabel('Time')
ylabel('Monomerization Fraction')
title("Monomerization Time Series with rbind="+rbind+" koff ="+koff)
ylim([0,1]);
axes('Position',[.7 .7 .2 .2])
box on
plot(tdat(1:5000),data(1:5000))
clear data;