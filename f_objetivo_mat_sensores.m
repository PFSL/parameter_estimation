function [fob] = f_objetivo_mat_sensores(par, nodes, elem, bc, acc, vel, des, alfa_damp, beta_damp, iddof, sdof, niddof, F, indexsensors, indexbc, E0, N, dt, ngl, nnel, nelem, nnos)

% placing young modulos to mount the system
for i=1:nelem
    elem(i,4) = E0*par(2+i);
end

% damping coefficients
alfa_d = par(1);
beta_d = par(2);

%% mount system using FEM
[Mef, Cef, Kef, Fef, iddof, indexnodof2, nidof, sdof, ngl] = truss_FEM(nodes, elem, bc, alfa_d, beta_d, N, F, nelem, ngl, sdof, nnos, nnel);

%% solve using state space equations
[y1,x1] = state_space_dyn(nidof, N, dt, Mef, Kef, Cef, Fef);

% output model
acc_inv = zeros(sdof,N);       % accelerations
vel_inv = zeros(sdof,N);       % velocities
des_inv = zeros(sdof,N);       % displacements

% indice vectors pointers
des_id = 1:nidof;
vel_id = nidof+1:2*nidof;
acc_id = 2*nidof+1:3*nidof;
Fob_id = 3*nidof+1:4*nidof;

% separated vectors
des_inv(iddof,:) = y1(des_id,:);
vel_inv(iddof,:) = y1(vel_id,:);
acc_inv(iddof,:) = y1(acc_id,:);

%% objective function ---> most important
% using only nodes with sensors
index_inv = indexsensors;

% objetive function based only in acceleration measures
fob = acc(index_inv,:) - acc_inv(index_inv,:);  % difference between accleerations

end