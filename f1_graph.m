%**************************************************************************

% ****** This function represent the two objective functions used here: one is
% contribution and another is centrality ******************************************************

%**************************************************************************

function f=f1_graph(c,adj,shortest_path_matrix)
% f1=(sum(sum(adj(c,c),1)))/((length(c))*(length(c))-(length(c)));

% f2=No_Of_Interaction_Partner(c,adj);    %
% f2_score=f2/length(c);
mat=adj(c,c);
v1=sum(mat,2);
v2=sum(adj(c,:),2);
contr_f2=sum(v1./v2);
% f=1/(f1+f2_score);
for i=1:length(c)
%     path_length=graphshortestpath(sparse(adj),C(i));
%     v=find(path_length==Inf);
%     path_length(v)=0;
%     centrality(i)=1/sum(path_length);

centrality(i)=1/sum(shortest_path_matrix(c(i),:));

end
f=contr_f2+ 1/sum(centrality);



