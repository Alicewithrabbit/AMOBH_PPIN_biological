%**************************************************************************

%This function performs selection, mutation etc. on a parent population and return the new child population ****************

%**************************************************************************
%**************************************************************************

function [New_Pop]=Make_New_Pop_PPIN_Biological(Pop_Size,Index_Of_Chromosome,adj,shortest_path_matrix,similarity_matrix)


[Sorted_Individual,F]=Non_Dom_Sort_PPIN_Biological(Pop_Size,Index_Of_Chromosome,adj,shortest_path_matrix,similarity_matrix);

for p=1:Pop_Size
    a=randperm(Pop_Size);
    c1=a(1);
    c2=a(2);

    if Sorted_Individual(c1).rank==Sorted_Individual(c2).rank

        [Sorted_Front_Ind,Distance_Of_Front]=Crowded_Dist_PPIN_Biological(F{Sorted_Individual(c1).rank},Index_Of_Chromosome,adj,shortest_path_matrix,similarity_matrix);
        for i=1:length(Sorted_Front_Ind)
            if c1==Sorted_Front_Ind(i)
                dist1=Distance_Of_Front(i);
            end
            if c2==Sorted_Front_Ind(i)
                dist2=Distance_Of_Front(i);
            end
        end
    end

    if (Sorted_Individual(c1).rank < Sorted_Individual(c2).rank) | ((Sorted_Individual(c1).rank==Sorted_Individual(c2).rank) && (dist1 > dist2))
        mat{p}=Index_Of_Chromosome{c1};

    else
        mat{p}=Index_Of_Chromosome{c2};
    end
end


New_Pop=[];

for i=1:Pop_Size
    a=randperm(Pop_Size);


    parent = mat{a(1)};

    if rand < 0.9

        [offspring]=change_chromosome1(parent,adj);
        New_Pop{i}=offspring;

    else
        New_Pop{i}=parent;


    end
end





