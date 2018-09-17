function [K] = asmbl(K, k_loc, elem, el, nnel, ngl)

%% mount assembly vector used to place local dof in global system
nodes_el = [elem(el,2), elem(el,3)];

% takes dofs of nodes in local system
k = 0;
for i=1:nnel
    ini = (nodes_el(i)-1)*ngl;

    for j=1:ngl
        k = k + 1;
        index(k) = ini + j;

    end
end

edof = length(index);

% converts in a position in global system 
for i=1:edof
    ii = index(i);
    
    for j=1:edof
        jj = index(j);
        K(ii,jj) = K(ii,jj) + k_loc(i,j);
        
    end
    
end
end