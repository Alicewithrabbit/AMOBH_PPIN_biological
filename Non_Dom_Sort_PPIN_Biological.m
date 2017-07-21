function [individual, Front]=Non_Dom_Sort_PPIN_Biological(Pop_Size,Index_Of_Chromosome,adj,shortest_path_matrix,similarity_matrix)
 
%pop0=round(rand(Pop_Size,8));
Front{1}=[];
for i=1:Pop_Size
    
    individual(i).Sp=[];
    individual(i).np=0;
    fv1=[f1_graph(Index_Of_Chromosome{i},adj,shortest_path_matrix),f2_Similarity(Index_Of_Chromosome{i},similarity_matrix)];%Ŀ�꺯��
    for j=1:Pop_Size
%         x1=dec(pop0(i,:));
%         x2=dec(pop0(j,:));
        %         if (obj1(x1)<obj1(x2) & obj2(x1)<=obj2(x2)) | (obj1(x1)<=obj1(x2) & obj2(x1)<obj2(x2))
         
        fv2=[f1_graph(Index_Of_Chromosome{j},adj,shortest_path_matrix),f2_Similarity(Index_Of_Chromosome{j},similarity_matrix)];
        
%         t= Check_Dominance_PPIN_Biological(fv1,fv2);     %individual(i) dominate individual(j) Check_Dominance(individual(i))
        if Check_Dominance_PPIN_Biological(fv1,fv2)  %���Paretoռ�ţ����i��ռ��j��
        individual(i).Sp=[individual(i).Sp,j];%��iռ�ŵĽ⼯
        else if  Check_Dominance_PPIN_Biological(fv2,fv1)         % ind(i)is dominated by ind(j)check_Dominance(%(obj1(x2)<obj1(x1) & obj2(x2)<=obj2(x1)) | (obj1(x2)<=obj1(x1) & obj2(x2)<obj2(x1))
 
                %                 Check_Dominance(scalling(dec(pop0(j,:))),scalling(dec(pop0(i,:))))
                individual(i).np=individual(i).np+1;%ռ��i��Ľ����
            end
        end
    end
% end
%k=1;
% for i=1:Pop_Size
 
    if individual(i).np==0%���������������1.i�������������⻥�಻ռ�ţ�2.i��ռ�����н⡣
       % indx(k)=i;
        individual(i).rank=1;
        Front{1}=[Front{1},i];%ǰ�˼���i��
    %k=k+1;
    end
end
 

%�� 
i=1;
while ~isempty(Front{i})
    Q=[];
    for j=1:length(Front{i})
        for k=1:length(individual(Front{i}(j)).Sp)
            individual(individual(Front{i}(j)).Sp(k)).np=individual(individual(Front{i}(j)).Sp(k)).np-1;
            if individual(individual(Front{i}(j)).Sp(k)).np==0
                individual(individual(Front{i}(j)).Sp(k)).rank=i+1;
                Q=[Q,individual(Front{i}(j)).Sp(k)];
            end
        end
    end
    i=i+1;
    Front{i}=Q;%�������ǰ�˷ּ�
end




