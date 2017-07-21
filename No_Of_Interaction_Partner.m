%*********This function compute the number of outward interaction partner
%for a chromosome c ******************************************************
%计算邻点数量
function [Interaction_Partner_Outward,r1]=No_Of_Interaction_Partner(c,adj)
I=[];
for i=1:length(c)
    s{i}=(find(adj(c(i),:)));
I=[I,s{i}];
end
Interaction_Partner_Total=unique(I);
common=wcommon(Interaction_Partner_Total,c);
indx=find(common==1);
[q,r]=qr(Interaction_Partner_Total);
[q1,r1]=qrdelete(q,r,indx);
Interaction_Partner_Outward=length(Interaction_Partner_Total)-length(c);