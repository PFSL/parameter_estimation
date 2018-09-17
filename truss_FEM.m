function [Mef, Cef, Kef, Fef, iddof, indexnodof2, nidof, sdof, ngl] = truss_FEM(nodes, elem, bc, alfa_damp, beta_damp, N, F, nelem, ngl, sdof, nnos, nnel)


%% allocate memory
K = zeros(sdof,sdof);       % stiffness
M = zeros(sdof,sdof);       % mass

%% mount global systema
[K, M] = global_sys(K, M, nodes, elem, nelem, ngl, nnel);

%% evaluate damping matriz using proportional rayleigh damping
C = alfa_damp.*M + beta_damp.*K;

%% effetive system
% apply bc's
[nbc, ~] = size(bc);

k = 0;
for i=1:nbc
    ini = (bc(i,1)-1)*ngl;          % find position at vector
    for j=1:ngl
        if bc(i,j+1) == 1
            k = k + 1;             
            indexbc(k) = ini+j;
        end
    end    
end

% select dof's to be evaluated 
k=0;
for i=1:nnos
    for j=1:ngl
        k=k+1;        
        indexnodof(k,1) = k;
        indexnodof(k,2) = nodes(i,1);
    end
end

indexnodof2 = indexnodof';

% vector of dof's to be solved
iddof = setdiff(indexnodof2(1,:),indexbc);
[~, nidof] = size(iddof);

% mount effective matrices to de used
Mef = M(iddof,iddof);
Kef = K(iddof,iddof);
Cef = C(iddof,iddof);
Fef = F(iddof,:);

end