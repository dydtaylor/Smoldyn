smoldynDat = {};
repoDat = {};
for i = 2:18
    repoFile = "SmoldynSimpleMFPT/Data/nSweep/"+i+".txt";
    repoDat{i-1} = load(repoFile);
end

for i = 1:17
    repoMeans(i) = mean(repoDat{i});
    repoErr(i)=std(repoDat{i})/sqrt(numel(repoDat{i}));
end
figure(2);
clf;
set(gca,'yscale','log')
hold on
errorbar(2:18,repoMeans,repoErr,'o','linewidth',2);
title('Mean First Passage Times')
xlabel('Number of Particles')
ylabel('MFPT')