% this function perform mutation on each chromosome.
% randomly select a node then performs either of the two task with equal
% probality 1.delete the node or
%           2.add the edges which are not included in this chromosome
%*********************************************************************

function [off]=change_chromosome1(parent,adj)
No_Of_Change=5;
for i=1:No_Of_Change
    a=randperm(length(parent));
    victim_node=parent(a(1));
    if rand < 0.5

        [q,r]=qr(parent);
        [q1,r1]=qrdelete(q,r,a(1));

        if length(r1)<=1

            m=find(adj(r1,:));
            off=union(r1,m);

        else

            off=r1;
        end

        parent=off;
    else

        m=find(adj(victim_node,:));
        off=union(parent,m);
        parent=off;
    end
end