%**************************************************************************
%********* This function formed the chromosome from biclustering result and
%calculates the density of each chrosome***********************************

%**************************************************************************

function d=f1_calc(c1,c2,adj)
c=union(c1,c2);
d=(sum(sum(adj(c,c),1)))/((length(c))*(length(c))-(length(c)));

 