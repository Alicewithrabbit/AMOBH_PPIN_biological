%**************************************************************************
%This function perform biclustering on a network specified by
%adjacency matrix (adj) and returns 'Pop_Size' number of highly dense
%bicluster after sorting each bicluster based on density. Here number of
%bicluster is specified as 'No_Of_Cluster*No_Of_Cluster'.

%**************************************************************************
%**************************************************************************

function [chrom]=biclustering(adj,No_Of_Cluster,Pop_Size)
fc1=kmeans(adj,No_Of_Cluster)';%kmeans���෵��������ǩ���������Ǳ������������ԡ�
fc2=kmeans(adj',No_Of_Cluster)';%kmeans���෵��������ǩ���������Ǳ������������ԣ�ת�ý���˫�����
l=1
for i=1:No_Of_Cluster
    for j=1:No_Of_Cluster
        a=find(fc1==i);
        b=find(fc2==j);
        if length(union(a,b))<3 %�����������۳���һ����ά���࣬�����ĸ���С����������
            continue;
        else
            SM(l).rows=a;
            SM(l).cols=b;

            dens(l)=f1_calc(a,b,adj);%�������ܶ�
            l=l+1;
        end
    end
end
[p,q]=sort(dens,'descend');
for i=1:Pop_Size
    chrom{i}= union(SM(q(i)).rows,SM(q(i)).cols);%���о������о���ȡ��������Ⱦɫ��
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
