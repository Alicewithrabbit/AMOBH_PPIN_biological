%********This function compute the crowded ditance between two individuals in a front*********

function [sorted_front,sorted_dist]=Crowded_Dist_PPIN_Biological(F,Index_Of_Chromosome,adj,shortest_path_matrix,similarity_matrix)
l=length(F);
% for i=1:l
%     dist1(i)=0;
%     dist2(i)=0;
% end

dist1 = zeros(1,l);
dist2 = zeros(1,l);

for i=1:l
    obj_value1(i)=f1_graph(Index_Of_Chromosome{F(i)},adj,shortest_path_matrix);%目标1
    obj_value2(i)=f2_Similarity(Index_Of_Chromosome{F(i)},similarity_matrix);%目标2
end

f1_max=max(obj_value1);
f1_min=min(obj_value1);
f2_max=max(obj_value2);
f2_min=min(obj_value2);
[x1,y1]=sort(obj_value1);
[x2,y2]=sort(obj_value2);
dist1(1)=99999;
dist1(l)=99999;
dist2(1)=99999;
dist2(l)=99999;
for i=2:(l-1)

    dist1(i)=dist1(i)+ abs((f1_graph(Index_Of_Chromosome{F(y1(i+1))},adj,shortest_path_matrix) - f1_graph(Index_Of_Chromosome{F(y1(i-1))},adj,shortest_path_matrix)))/(f1_max-f1_min);
    dist2(i)=dist2(i)+abs((f2_Similarity(Index_Of_Chromosome{F(y2(i+1))},similarity_matrix) - f2_Similarity(Index_Of_Chromosome{F(y2(i-1))},similarity_matrix)))/(f2_max-f2_min);
end

for i=1:l
    dist(i)=dist1(i)+dist2(i);
end
[sorted_dist,indx]=sort(dist,2,'descend');
for i=1:l
    sorted_front(i)=F(indx(i));
end