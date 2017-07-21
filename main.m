protein_A=MOESMESM(:,1);
protein_B=MOESMESM(:,2);
Humain_proteins_list = unique(union(protein_A, protein_B));
Humain_proteins_data = sparse(length(Humain_proteins_list),length(Humain_proteins_list));
for i = 1:length(protein_A)
     i_indices = strcmp(Humain_proteins_list,protein_A(i));
     j_indices = strcmp(Humain_proteins_list,protein_B(i));
     Humain_proteins_data(i_indices,j_indices) = Humain_proteins_data(i_indices,j_indices)+1;
     Humain_proteins_data(j_indices,i_indices) = Humain_proteins_data(i_indices,j_indices);
 end
Humain_proteins_data(Humain_proteins_data>0)=1;
no_weight_ppi_data = Humain_proteins_data;
no_weight_ppi_data(no_weight_ppi_data>0) = 1;
[S,C]=graphconncomp(logical(sparse(no_weight_ppi_data)),'Directed','false');
connected_proteins_list_index = find(C==mode(C));
Humain_connected_proteins_data = Humain_proteins_data(connected_proteins_list_index,connected_proteins_list_index);
Humain_connected_ppi_proteins_list = Humain_proteins_list(connected_proteins_list_index);
save (['C:\Users\Administrator\Desktop\amobh_protein\','Humanpxy','_data'],'Humain_proteins_data','Humain_proteins_list','Humain_connected_proteins_data', 'Humain_connected_ppi_proteins_list')
%将数据存储起来