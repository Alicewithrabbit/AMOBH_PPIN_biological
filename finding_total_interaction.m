function l=finding_total_interaction(complex_protein,interaction_data) 
 for i=1:length(complex_protein)
     
[a1,b1]=intersect(complex_protein{i},interaction_data(:,1));
l1=length(intersect(complex_protein{i},interaction_data(b1,2)));
[a2,b2]=intersect(complex_protein{i},interaction_data(:,2));
l2=length(intersect(complex_protein{i},interaction_data(b2,1)));
l(i)=l1+l2;
 end