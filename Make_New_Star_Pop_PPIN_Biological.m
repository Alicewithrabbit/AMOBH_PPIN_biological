function [New_Pop]= Make_New_Star_Pop_PPIN_Biological(Index_Of_Chromosome,Old_Pop,adj,nvars)

Parent_Pop = Index_Of_Chromosome;
New_Pop = [];

for i = 1:nvars
    
    if rand < 0.9
    a=randperm(size(Parent_Pop,2));
    
    selected_bh = a(1);
    
    black_hole = Parent_Pop{selected_bh};
    

    [offspring]=change_chromosome1(black_hole,adj);
    
    New_Pop{i}=offspring;
    
    else
        
        a=randperm(size(Old_Pop,2));
        
        
        [offspring]=change_chromosome1(Old_Pop{a(1)},adj);
    
        New_Pop{i}=offspring;
        
        
    end
end

end





