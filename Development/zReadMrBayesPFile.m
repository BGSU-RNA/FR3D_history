
% t= textread('23s_data_19bins.nex.run2.p','','headerlines',2);
t= textread('5s_data_19bins.nex.run1.p','','headerlines',2);

m = mean(t((end-699):end,10:end));

save 5srates.mat m
