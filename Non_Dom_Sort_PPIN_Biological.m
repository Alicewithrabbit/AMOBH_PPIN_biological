function [individual, Front]=Non_Dom_Sort_PPIN_Biological(Pop_Size,Index_Of_Chromosome,adj,shortest_path_matrix,similarity_matrix)
 
%pop0=round(rand(Pop_Size,8));
Front{1}=[];
for i=1:Pop_Size
    
    individual(i).Sp=[];
    individual(i).np=0;
    fv1=[f1_graph(Index_Of_Chromosome{i},adj,shortest_path_matrix),f2_Similarity(Index_Of_Chromosome{i},similarity_matrix)];%目标函数
    for j=1:Pop_Size
%         x1=dec(pop0(i,:));
%         x2=dec(pop0(j,:));
        %         if (obj1(x1)<obj1(x2) & obj2(x1)<=obj2(x2)) | (obj1(x1)<=obj1(x2) & obj2(x1)<obj2(x2))
         
        fv2=[f1_graph(Index_Of_Chromosome{j},adj,shortest_path_matrix),f2_Similarity(Index_Of_Chromosome{j},similarity_matrix)];
        
%         t= Check_Dominance_PPIN_Biological(fv1,fv2);     %individual(i) dominate individual(j) Check_Dominance(individual(i))
        if Check_Dominance_PPIN_Biological(fv1,fv2)  %检查Pareto占优，如果i解占优j解
        individual(i).Sp=[individual(i).Sp,j];%被i占优的解集
        else if  Check_Dominance_PPIN_Biological(fv2,fv1)         % ind(i)is dominated by ind(j)check_Dominance(%(obj1(x2)<obj1(x1) & obj2(x2)<=obj2(x1)) | (obj1(x2)<=obj1(x1) & obj2(x2)<obj2(x1))
 
                %                 Check_Dominance(scalling(dec(pop0(j,:))),scalling(dec(pop0(i,:))))
                individual(i).np=individual(i).np+1;%占优i解的解个数
            end
        end
    end
% end
%k=1;
% for i=1:Pop_Size
 
    if individual(i).np==0%这里有两种情况，1.i解与所有其他解互相不占优，2.i解占优所有解。
       % indx(k)=i;
        individual(i).rank=1;
        Front{1}=[Front{1},i];%前端加入i解
    %k=k+1;
    end
end
 

%对 
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
    Front{i}=Q;%按排序给前端分级
end




