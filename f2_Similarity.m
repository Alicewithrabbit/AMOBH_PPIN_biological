%***********************************************************************
%This function represent 2nd objective used here. It compute the avg. semantic similarity of a chromosome c using similarity Matrix
%***********************************************************************

function obj=f2_Similarity(C,similarity_matrix)
f=mean(mean(similarity_matrix(C,C)));    
obj=1/f;

