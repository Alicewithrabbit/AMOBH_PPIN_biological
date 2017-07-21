%This is the main function. It performs a multi-objective GA based
%clustering in Protein interaction network.

%Input: Population size, No. of Generation, Adjacency Matrix of the whole
%network, 
%Similarity Matrix:containing the GO based semantic similarity between each
%pair of nodes
%shortest_path_matrix: containing the shortest path distance between each
%pair of node in the whole PPI network

% Output: New_Pop is the resulting population which consists of predicted
% protein complexes. It consists of 50 (Pop_Size) number of cell array. Each cell array consists of set of index of proteins.

%**************************************************************************
%**************************************************************************

function [front,New_Pop]=MOEA_PPIN_Biological(adj,Pop_Size,No_Of_Generation,shortest_path_matrix,similarity_matrix)

%**** Perform biclustering on the whole network and picking up 50 (Pop_Size) number of bicluster based on density***********************

chrom = biclustering(adj,Pop_Size/2,Pop_Size);
% load chrom.mat
%**************************************************************************

%****************** Making Initial Population *****************************

Initial_Pop = chrom;

pop0 = Make_New_Pop_PPIN_Biological(Pop_Size,Initial_Pop,adj,shortest_path_matrix,similarity_matrix)

%**************************************************************************

%********perform multiobjective optimization based on NSGA-II************** 
%%
for t=1:No_Of_Generation
    display('generation;');
    t
    disp('is created');

    Parent_Pop={};


    R=[Initial_Pop,pop0];

    [Sort_ind,F]=Non_Dom_Sort_PPIN_Biological(size(R,2),R,adj,shortest_path_matrix,similarity_matrix);

    i=1;
    while (size(Parent_Pop,2)+length(F{i})) <= Pop_Size



        for j=1:length(F{i})

            Parent_Pop=[Parent_Pop,R{F{i}(j)}];
        end
        i=i+1;
    end
    size(Parent_Pop,2)
    sorted_front=Crowded_Dist_PPIN_Biological(F{i},R,adj,shortest_path_matrix,similarity_matrix);
    for k=1:(Pop_Size - size(Parent_Pop,2))

        Parent_Pop=[Parent_Pop,R{sorted_front(k)}];
    end

    New_Pop=Make_New_Pop_PPIN_Biological(Pop_Size,Parent_Pop,adj,shortest_path_matrix,similarity_matrix);


    pop0=Parent_Pop;
    Initial_Pop=New_Pop;


end
[individual,front]=Non_Dom_Sort_PPIN_Biological(Pop_Size,New_Pop,adj,shortest_path_matrix,similarity_matrix);


%%
%**************************************************************************

%***********Drawing the front**********************************************

for i=1:length(front{1})
    solx(i)=f1_graph(New_Pop{front{1}(i)},adj,shortest_path_matrix);
    soly(i)=f2_Similarity(New_Pop{front{1}(i)},similarity_matrix);
end
plot(solx,soly,'*');

%**************************************************************************
