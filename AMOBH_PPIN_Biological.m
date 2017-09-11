%%
%AMOBH main function
%Author: Chong Wu
%Date: 2017.07.21
%Syntax: [front, gArchive] = AMOBH_PPIN_Biological(adj, shortest_path_matrix, similarity_matrix, bh_option)
%bh_option is a stucture.
%It performs an adaptive multi-objective black hole algorithm based
%on clustering in Protein interaction network.
function [front, gArchive] = AMOBH_PPIN_Biological(adj, shortest_path_matrix, similarity_matrix, bh_option)

Pop_Size = bh_option.sizestar;

Archive_Size = bh_option.narc;

No_Of_Generation = bh_option.maxgen;

ELS_MAX = bh_option.els_max;%Maximum of learning rate;

ELS_MIN = bh_option.els_min;%Minimum of learning rate;

ELS = ELS_MAX;

step_x = (ELS_MAX - ELS_MIN)/bh_option.maxgen;%step size of learning;

delta_s = 2/(Archive_Size*Pop_Size)*log2(2);

status = 'convergence';

gArchive = [];

nobjs = 2;

KT = 2*(nobjs+1);

%**** Perform biclustering on the whole network and picking up 50 (Pop_Size) number of bicluster based on density***********************
% 
% chrom = biclustering(adj,Pop_Size/2,Pop_Size);

%**************************************************************************
%****************** Making Initial Population *****************************
load chrom.mat
Initial_Pop = chrom;
[initial_archive, Entropy, delta_Entropy, gbest_all,L] = Initial_Archive(Pop_Size,Initial_Pop,adj,shortest_path_matrix,similarity_matrix, bh_option);
gArchive = initial_archive;

Entropy_temp = Entropy(1);

K = Archive_Size;

n1 = size(gArchive,2);

New_Pop = Initial_Pop;

%%
%迭代寻优
for t=1:No_Of_Generation
    
    display('generation;');
    t
    disp('is created');

    ngb = size(gbest_all,2);
    
    selector = randperm(ngb);
    
%     black_hole = gbest_all{selector};
    
    if rand < ELS
        
        gbest_all{selector(1)} = change_chromosome1(gbest_all{selector(1)},adj);
%         black_hole = gbest_all{selector(1)};
        
    end
    
    
    New_Pop = Make_New_Star_Pop_PPIN_Biological(gbest_all,New_Pop,adj,Pop_Size);
    

    
    fv1 = [];
    fv2 = [];
    
    
    
    
    diff = [];
    for i = 1:Pop_Size
        
        fv1(i,:)=[f1_graph(New_Pop{i},adj,shortest_path_matrix),f2_Similarity(New_Pop{i},similarity_matrix)];
        
        ngArc = size(gArchive,2);
        count = 0;
        I = 0;
        if ngArc == 1
            j = 1;
            fv2(j,:) = [f1_graph(gArchive{j},adj,shortest_path_matrix),f2_Similarity(gArchive{j},similarity_matrix)];
            if sum(all(fv1(i,:) <= fv2(j,:) ,1))>0 &&sum(any(fv1(i,:) < fv2(j,:),1))>0&&count == 0
                gArchive{j} = New_Pop{i};
                count = 1;

            end
            
        else
            
            for j = 1:ngArc
                fv2(j,:) = [f1_graph(gArchive{j},adj,shortest_path_matrix),f2_Similarity(gArchive{j},similarity_matrix)];
                if all(fv1(i,:) <= fv2(j,:)) &&any(fv1(i,:) < fv2(j,:))&&count == 0
                    gArchive{j} = New_Pop{i};
                    count = 1;
                    I(j) = 0;
                elseif  all(fv1(i,:) <= fv2(j,:))&&any(fv1(i,:) < fv2(j,:))&&count == 1


                    I(j) = 1;
                end
            end
        end
        
        p = find(I>0);
        if ~isempty(p)
            gArchive(p) =[];
            ngArc = size(gArchive,2);
        end
        flag = 0;
        fv2 = [];
        
        
        if count == 0
            if ngArc < Archive_Size
                
                if ngArc == 1
                    j = 1;
                        fv2(j,:) = [f1_graph(gArchive{j},adj,shortest_path_matrix),f2_Similarity(gArchive{j},similarity_matrix)];
                        any((fv1(i,:) <= fv2(j,:)))
                        if any((fv1(i,:) <= fv2(j,:)))
                              gArchive{ngArc + 1} = New_Pop{i};

                        end                    
                else
                    for j = 1:ngArc
                        fv2(j,:) = [f1_graph(gArchive{j},adj,shortest_path_matrix),f2_Similarity(gArchive{j},similarity_matrix)];
                        if any((fv1(i,:) <= fv2(j,:)))

                        elseif all((fv1(i,:) > fv2(j,:)))

                            flag = 1;

                        end
                    end
                if flag == 0
                    gArchive{ngArc + 1} = New_Pop{i};
                end                    
                    
                end
                

                ngArc = size(gArchive,2);
                
            else
                
                L_temp = [];
                gArchive_temp = [];
                gArchive_temp = [gArchive,New_Pop{i}];
                fv2 = [];
                
                for j = 1:ngArc + 1
                    fv2(j,:) = [f1_graph(gArchive{j},adj,shortest_path_matrix),f2_Similarity(gArchive{j},similarity_matrix)];
                end
                
                for j = 1:ngArc+1
                    
                    L_temp(j,1) = ceil(n*(fv2(j,1) - min(fv2(:,1)))/(max(fv2(:,1)) -...
                        min(fv2(:,1))));
                    L_temp(j,2) = ceil(n*(fv2(j,2) - min(fv2(:,2)))/(max(fv2(:,2)) -...
                        min(fv2(:,2))));
                    if fv2(j,1)  == min(fv2(:,1) )
                        
                        L_temp(j,1) = 1;
                        
                    end
                    if fv2(j,2)  == min(fv2(:,2) )
                        
                        L_temp(j,2) = 1;
                        
                    end
                end
                for j = 1:ngArc + 1
                    
                    den(j) = sum(sum(repmat(L_temp(j,:),size(L_temp,1),1) == L_temp));
                    
                end
                
                
                c = find(den == max(den));
                %若密度小于最大密度，进入文档
                select = c(randperm(length(c),1));
                if den(end) < den(select)
                    gArchive{select} = New_Pop{i};
                end
                
                
                
            end
            
        end
        ngArc = size(gArchive,2);
        
        
    end
    
    ngArc = size(gArchive,2);
    fv3 = [];
    L = [];
    for j = 1:ngArc
        fv3(j,:) = [f1_graph(gArchive{j},adj,shortest_path_matrix),f2_Similarity(gArchive{j},similarity_matrix)];
    end
    
    for j = 1:ngArc
        
        L(j,1) = ceil(ngArc*(fv3(j,1) - min(fv3(:,1)))/(max(fv3(:,1)) -...
            min(fv3(:,1))));
        L(j,2) = ceil(ngArc*(fv3(j,2) - min(fv3(:,2)))/(max(fv3(:,2)) -...
            min(fv3(:,2))));
        if fv3(j,1)  == min(fv3(:,1) )
            
            L(j,1) = 1;
            
        end
        if fv3(j,2)  == min(fv3(:,2) )
            
            L(j,2) = 1;
            
        end
    end
    
    for i = 1:nobjs
        for j = 1:ngArc
            if(~isempty(find(L(:,i) == j, 1)))
                Entropy(t) = Entropy(t) - length(find(L(:,i) == j))/ngArc*nobjs*log2(length(find(L(:,i) == j))/(ngArc*nobjs));
            else
                Entropy(t) = Entropy(t) - 0;
            end
        end
    end
    
    delta_Entropy(t) = Entropy(t) - Entropy_temp;
    
    n = ngArc;
    
    % check the status
    delta_c = 2/n*log2(2);
    Entropy_temp = Entropy(t);
    if abs(delta_Entropy(t)) > delta_c || abs(n - n1) >0
        status = 'convergence';
    end
    
    if abs(delta_Entropy(t)) < delta_c && n == K && n1 == K
        if abs(delta_Entropy(t)) > delta_s
            status = 'diversity';
        end
    end
    
    if abs(delta_Entropy(t)) < delta_s && n == n1
        status = 'stagnation';
    end
    
    
    den1 = [];
    if n~=1
        for i = 1:n
            
            den1(i) = sum(sum(repmat(L(i,:),size(L,1),1) == L));
            
        end
    else
        
        den1 = 3;
        
    end
    
    %状态收敛时
    if(strcmp(status,'convergence'))
        
        
        sc = [];
        scc = [];
        den2 = [];
        gbest_all = [];
        n1 = n;
        for k = 1:n
            sc(k) = score(L,n,k);
        end
        if length(den1) < KT
            
            gbest_all = gArchive;
            
        else
            den2 = sort(den1);
            for i = 1:KT/2-1
                I = find(den1 == den2(i),1);
                gbest_all{i} =  gArchive{I};
            end
            scc = sort(sc,'descend');
            for k = 1:KT/2+1
                I = find(sc == scc(k),1);
                gbest_all{k+KT/2-1} =  gArchive{I};
            end
            
        end
        
        
        
    end
    %状态多样时
    if(strcmp(status,'diversity'))
        ELS = ELS - step_x*abs(delta_Entropy(t));
        sc = [];
        scc = [];
        den2 = [];
        gbest_all = [];
        n1 = n;
        for k = 1:n
            sc(k) = score(L,n,k);
        end
        if length(den1) < KT
            
            gbest_all = gArchive;
            
        else
            den2 = sort(den1);
            for i = 1:KT/2+1
                I = find(den1 == den2(i),1);
                gbest_all{i} =  gArchive{I};
            end
            scc = sort(sc,'descend');
            for k = 1:KT/2-1
                I = find(sc == scc(k),1);
                gbest_all{k+KT/2+1} =  gArchive{I};
            end
            
        end
        
        
        
    end
    
    %状态停滞时
    if(strcmp(status,'stagnation'))
        ELS = ELS + 2*step_x*abs(1 + delta_Entropy(t));
        sc = [];
        scc = [];
        den2 = [];
        gbest_all = [];
        n1 = n;
        for k = 1:n
            sc(k) = score(L,n,k);
        end
        if length(den1) < KT
            
            gbest_all = gArchive;
            
        else
            den2 = sort(den1);
            for i = 1:KT/2
                I = find(den1 == den2(i),1);
                gbest_all{i} =  gArchive{I};
            end
            scc = sort(sc,'descend');
            for k = 1:KT/2
                I = find(sc == scc(k),1);
                gbest_all{k+KT/2} =  gArchive{I};
            end
            
        end
        
        
    end
    
    fv3 = [];
    
    for i = 1:size(gbest_all,2)
        fv3(i,:) = [f1_graph(gbest_all{i},adj,shortest_path_matrix),f2_Similarity(gbest_all{i},similarity_matrix)];
    end
    
    v = size(gbest_all,2);
    R = [];
    
    %计算黑洞半径
    if v == 1
        R = abs(fv3./(sum(fv1)));
    else
        for j = 1:v
            R = [R;abs(fv3(j)./sum(fv1))];
        end
    end
    
    % 判断是否进入黑洞边界
    for j = 1:bh_option.sizestar
        %         S = min(abs(solution_vector(j,:) - gbest(1,1:M)));
        if any(sum(abs(repmat(fv1(j,:),size(gbest_all,2),1) - fv3) > R) ==0) %进入吸收边界:吸收并重新生成
            
            sel = randperm(size(chrom,2));
            New_Pop{j} = chrom{sel(1)};
            
        end
    end
    
    
end

%***********Drawing the front**********************************************



for i=1:size(gArchive,2)
    solx(i)=f1_graph(gArchive{i},adj,shortest_path_matrix);
    soly(i)=f2_Similarity(gArchive{i},similarity_matrix);
end
front = [solx;soly];
plot(solx,soly,'*');
