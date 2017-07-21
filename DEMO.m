%%
%AMOBH 
% bh_option = struct('maxgen',60,'sizestar',50,'els_max',0.2,'els_min',0.1,'narc',50);
% [front,gArchive] = AMOBH_PPIN_Biological(adj,shortestpath,SIM,bh_option)
%%
%NSGA-II
[front,New_Pop]=MOEA_PPIN_Biological(adj,50,10,shortestpath,SIM)