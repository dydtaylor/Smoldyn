function mfptECDF(n)
file = "SmoldynSimpleMFPT/Data/nSweep/"+n+".txt";
data = load(file);
figure()
[f,x] = ecdf(data);
plot(x,f)
xlabel('Time')
ylabel('Cumulative Fraction')
titlestr = "Empirical CDF of First Passage Times for " + n + " Particles";
title(titlestr)