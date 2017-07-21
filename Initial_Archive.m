function [initial_archive, Entropy, delta_Entropy, gbest_all, L] = Initial_Archive(Pop_Size,Index_Of_Chromosome,adj,shortest_path_matrix,similarity_matrix, bh_option)

[Sorted_Individual,F,fv] = NDS_PPIN_Biological(Pop_Size,Index_Of_Chromosome,adj,shortest_path_matrix,similarity_matrix);


initial_archive_size = size(F{1}, 2);

if initial_archive_size < bh_option.narc

    for i = 1:initial_archive_size

        initial_archive{i} = Index_Of_Chromosome{F{1}(i)};

    end

else
    
    for i = 1:bh_option.narc

        initial_archive{i} = Index_Of_Chromosome{F{1}(i)};

    end
    
end

%%PCCS
nobjs = 2;
for i = 1:initial_archive_size

	for j = 1:nobjs

		L(i,j) = ceil(initial_archive_size*(fv(F{1}(i),j) - min(fv(F{1}(:),j)))/(max(fv(F{1}(:),j)) -...
			min(fv(F{1}(:),j))));
	
		if fv(F{1}(i),j) == min(fv(F{1}(:),j))

			L(i,j) = 1;
		
		end

	end

end

%%计算香浓信息熵
Entropy  = zeros(1,bh_option.maxgen);

Entropy_temp = Entropy(1);

for t = 1:nobjs
    for k = 1:initial_archive_size
        if(~isempty(find(L(:,t) == k, 1)))
            Entropy(1) = Entropy(1) - length(find(L(:,t) == k))/initial_archive_size*nobjs*log2(length(find(L(:,t) == k))/(initial_archive_size*nobjs));
        else
            Entropy(1) = Entropy(1) - 0;
        end
    end
end

delta_Entropy = Entropy - Entropy_temp;




%%找出初始黑洞集
sc = [];
gbest_all = [];
n1 = initial_archive_size;
if n1 == 1
    gbest_all = initial_archive;
else
    for k = 1:n1
        sc(k) = score(L,n1,k);
    end
    if all(sc == 0)
        gbest_all = initial_archive;
    else
        q = find(sc == max(sc));
%         sel = randperm(length(q));
        if length(q) == 1
            gbest_all = initial_archive{q(1)};
        else
            
            for j = 1:length(q)
                gbest_all{j} =  initial_archive{q(1)};
            end
        end
    end        
end
