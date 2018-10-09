function [K, M] = global_sys(K, M, nodes, elem, nelem, ngl, nnel)

%% mounting global stiffness and mass matrix using Finite Element Method

for el=1:nelem
    % physical properties
    rho = elem(el,6);
    A = elem(el,5);
    E = elem(el,4);
    
    % nodal coordinates of element
    x1 = nodes(elem(el,2),2);
    x2 = nodes(elem(el,3),2);
    y1 = nodes(elem(el,2),3);
    y2 = nodes(elem(el,3),3);
    
    % element length
    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    
    % element angle with global system
    if (x2-x1)==0
        beta = 2*atan(1);
    else
        beta = atan((y2-y1)/(x2-x1));
    end
        
    % local mass and stiffness evaluation
%     ipt = 1; % consistent
    ipt = 2; % lumped
    [k_loc, m_loc] = stiff_local2d(L, A, E, rho, beta, ipt);
    
    % assembly matrix
    [K] = asmbl(K, k_loc, elem, el, nnel,ngl);
    [M] = asmbl(M, m_loc, elem, el, nnel,ngl);
    
end
end