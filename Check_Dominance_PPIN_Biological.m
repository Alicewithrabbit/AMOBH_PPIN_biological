%**************************************************************************
% This function checks the dominating the dominating properties of two
% individuals
%**************************************************************************
%**************************************************************************

function t=Check_Dominance_PPIN_Biological(fv1,fv2)

t=0;%

no=size(fv1,2);

diff=fv1(1,:)-fv2(1,:);
s1=size(find(diff>0),2);
s2=size(find(diff<0),2);

if s2==no|(s2>0 & s1==0) %1占优2
   t=1;
end

if s1>0 & s2>0 % 互相不占优或2占优1
   t=0;
end

