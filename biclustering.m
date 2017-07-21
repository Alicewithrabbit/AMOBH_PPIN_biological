%**************************************************************************
%This function perform biclustering on a network specified by
%adjacency matrix (adj) and returns 'Pop_Size' number of highly dense
%bicluster after sorting each bicluster based on density. Here number of
%bicluster is specified as 'No_Of_Cluster*No_Of_Cluster'.

%**************************************************************************
%**************************************************************************

function [chrom]=biclustering(adj,No_Of_Cluster,Pop_Size)
fc1=kmeans(adj,No_Of_Cluster)';%kmeans聚类返回隶属标签向量，行是变量，列是属性。
fc2=kmeans(adj',No_Of_Cluster)';%kmeans聚类返回隶属标签向量，行是变量，列是属性，转置进行双向聚类
l=1
for i=1:No_Of_Cluster
    for j=1:No_Of_Cluster
        a=find(fc1==i);
        b=find(fc2==j);
        if length(union(a,b))<3 %两个聚类所聚出的一个二维聚类，该类点的个数小于三，舍弃
            continue;
        else
            SM(l).rows=a;
            SM(l).cols=b;

            dens(l)=f1_calc(a,b,adj);%排序用密度
            l=l+1;
        end
    end
end
[p,q]=sort(dens,'descend');
for i=1:Pop_Size
    chrom{i}= union(SM(q(i)).rows,SM(q(i)).cols);%对行聚类与列聚类取并集构成染色体
end



% function [chrom]=biclustering(adj,No_Of_Cluster,Pop_Size)
% fc1=kmeans(adj,No_Of_Cluster)';
% fc2=kmeans(adj',No_Of_Cluster)';
% for i=1:No_Of_Cluster
%     for j=1:No_Of_Cluster
%         a=find(fc1==i);
%         b=find(fc2==j);
%         if length(union(a,b))<3
%             continue;
%         else
%             SM(l).rows=a;
%             SM(l).cols=b;
% 
%             dens(l)=f1_calc(a,b,adj);
%             l=l+1;
%         end
%     end
% end
% [p,q]=sort(dens,'descend');
% for i=1:Pop_Size
%     chrom{i}=union(SM(q(i)).rows,SM(q(i)).cols);
% end
