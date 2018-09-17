clear all; clc; %close all

%% input of a simple truss
[nodes, elem, bc] = mimo();
% [nodes, elem, bc] = input_truss_0();

% define size of the problem
[nnos,~] = size(nodes); % number of nodes
[nelem,~] = size(elem); % number of elements
ngl = 2;                % number of degrees of freedom (dof) per node
sdof = nnos*ngl;        % total number of dof
nnel = 2;               % number of nodes per element

% define time properties
Fs = 15;            % Sample frequency (Hz)
dt = 1/Fs;          % time intervl (sec)
ti = 0;             % initial time
tf = 30;           % final time

t = [ti:dt:tf+dt]';    % time vector (tf+dt in order to consider initial time)
N = size(t,1);     % step times (minos 1 in order to not count initial condition)

% %% input truss design check
% plot_graph(nodes, elem, ngl)

%% load vector (varying in time)
F = zeros(sdof,N);

% for the specifc case
% F(5,:) = 10;
fhz1 = .10;                         % frequency 1 in Hz
fhz2 = .16;                         % frequency 2 in Hz
P0 = [10.0, 8.0];					% force amplitude
W0 = [fhz1*2*pi, fhz2*2*pi];		% frequency of excitation in rad/sec

for k = 1:N
    F(3,k) = P0(1)*sin(W0(1)*(k)*dt);
%     F(3,k) = P0(1)*sin(W0(1)*(k)*dt) + P0(2)*cos(W0(2)*(k)*dt);
    F(5,k) = P0(2)*sin(W0(2)*(k)*dt);
    
end

% figure()
% stem(t, F(3,:), 'filled', 'MarkerSize', 4)

%% proportional damping coefficients
alfa_damp = .10;
beta_damp = .10;

%% mounting the system trought FEM
[Mef, Cef, Kef, Fef, iddof, indexnodof2, nidof, sdof, ngl] = truss_FEM(nodes, elem, bc, alfa_damp, beta_damp, N, F, nelem, ngl, sdof, nnos, nnel);

%% numerical solver in state space
[y,x] = state_space_dyn(nidof, N, dt, Mef, Kef, Cef, Fef);

%% mount independent vectors from observer
acc = zeros(sdof,N);       % acelerations
vel = zeros(sdof,N);       % velocities
des = zeros(sdof,N);       % displacements
Fob = zeros(sdof,N);       % forces

% indice vectors pointers
des_id = 1:nidof;
vel_id = nidof+1:2*nidof;
acc_id = 2*nidof+1:3*nidof;
Fob_id = 3*nidof+1:4*nidof;

% separated vectors
des(iddof,:) = y(des_id,:);
vel(iddof,:) = y(vel_id,:);
acc(iddof,:) = y(acc_id,:);
Fob(iddof,:) = y(Fob_id,:);


%% plot results from system with 2 mass (mimo example)
figure()
% force in mass 1
subplot(4,2,1)
plot(t,F(3,:))
xlim([0 50])
ylim([-15 15])
title('Force at mass 1')
xlabel('time')
ylabel('F')

% force in mass 2
subplot(4,2,2)
plot(t,F(5,:))
xlim([0 50])
ylim([-15 15])
title('Force at mass 2')
xlabel('time')
ylabel('F')

% displacement in mass 1
subplot(4,2,3)
plot(t,des(3,:),'b')
xlim([0 50])
ylim([-50 50])
title('Displacement at mass 1')
xlabel('time')
ylabel('u')

% displacement in mass 2
subplot(4,2,4)
plot(t,des(5,:),'b')
xlim([0 50])
ylim([-50 50])
title('Displacement at mass 1')
xlabel('time')
ylabel('u')

% velocity in mass 1
subplot(4,2,5)
plot(t,vel(3,:),'b')
xlim([0 50])
ylim([-50 50])
title('Velocity at massa 1')
xlabel('time')
ylabel('v')

% velocity in mass 2
subplot(4,2,6)
plot(t,vel(5,:),'b')
xlim([0 50])
ylim([-50 50])
title('Velocity at massa 1')
xlabel('time')
ylabel('v')

% acceleration in mass 1
subplot(4,2,7)
plot(t,acc(3,:),'b')
xlim([0 50])
ylim([-50 50])
title('Acceleration at massa 1')
xlabel('time')
ylabel('a')

% acceleration in mass 2
subplot(4,2,8)
plot(t,acc(5,:),'b')
xlim([0 50])
ylim([-50 50])
title('Acceleration at massa 1')
xlabel('time')
ylabel('a')

%% inversion to obtain physical parameters of system given a output response
% mont effective system (removes bc's)
[nbc, ~] = size(bc);
k = 0;
n = 0;
for i=1:nbc
    ini = (bc(i,1)-1)*ngl;          % find position at vetor
    for j=1:ngl
        if bc(i,j+1) == 1
            k = k + 1;
            indexbc(k) = ini+j;
        end
    end    
end

% select dof's to be evaluated during inversion
k=0;
for i=1:nnos
    for j=1:ngl
        k=k+1;
        indexnodof(k,1) = k;
        indexnodof(k,2) = nodes(i,1);
    end
end


indexnodof2 = indexnodof';
% removes from all dofs in a vector those from bc's
iddof = setdiff(indexnodof2(1,:),indexbc);

[~, niddof] = size(iddof);

% sensors choosing (1-measured, 0-ignored)
% [sensor number, node measured, x direction, y direction]
sensors = [
%     1, 1, 1, 1;
    2, 2, 1, 0;
    3, 3, 1, 0;
%     4, 4, 1, 1;
    ];

% % plot nodes to be measured, sensors position and bc's
% plot_check_inv_prob(nodes, elem, ngl, iddof, bc, sensors);

[nsensor, ~] = size(sensors);
[nbc, lixo] = size(bc);

% select dof's measureds by sensors to be used in f_objetivo
k=0;
for i=1:nsensor
    ini = (sensors(i,2)-1)*ngl;          % find position at vector
    for j=1:ngl
        if sensors(i,j+2) == 1
            k = k + 1;
            indexsensors(k) = ini+j;
        end
    end    
end


%% identifica material das barras sem amortecimento
par0 = ones(nelem+2,1);     % Young parameters and damping coefficients
E0 = 5e0;                   % initial Young modulus
par0(:,1) = par0(:,1)*.1;   % multiplier used for estimate young modulus

alfa0 = 1.001;              % initial alpha damping 
beta0 = 1.001;              % initial beta damping
par0(1,1) = alfa0;          % set at par0 alfa0
par0(2,1) = beta0;          % set at par0 beta0

%% define objective function to estimate the inverse solution
f = @(x)f_objetivo_mat_sensores(x, nodes, elem, bc, acc, vel, des, alfa_damp, beta_damp, iddof, sdof, niddof, F, indexsensors, indexbc, E0, N, dt, ngl, nnel, nelem, nnos)

lb = par0*(0.0001); % lower boundary
ub = par0*(5.);     % upper boundary

% options used for lsqnonlin function
options = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective', 'MaxFunctionEvaluations',2000,'MaxIterations',12000, 'StepTolerance',1e-12, 'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12, 'PlotFcn',{@optimplotstepsize,@optimplotx,@optimplotresnorm}, 'Display', 'final-detailed', 'FiniteDifferenceType', 'forward');

%% call inversion routine as an optimization problem
[x,fval] = lsqnonlin(f,par0,lb,ub,options)

% show results
elem_inv = elem;
[nelem, ~] = size(elem);

% set parameters obtained in element matrix
for i=1:nelem
    elem_inv(i,4) = E0*x(2+i);
end

% show element matrix and dampings obtained 
elem_inv
alfa_damp_inv = x(1)
beta_damp_inv = x(2)


%% generate solution with parameters estimated and check with the original ones
[Mef, Cef, Kef, Fef, iddof, indexnodof2, nidof, sdof, ngl] = truss_FEM(nodes, elem, bc, alfa_damp, beta_damp, N, F, nelem, ngl, sdof, nnos, nnel);

[yinv,xinv] = state_space_dyn(nidof, N, dt, Mef, Kef, Cef, Fef);

acc_invertido = zeros(sdof,N);       % aceleration vector
vel_invertido = zeros(sdof,N);       % velocities vector
des_invertido = zeros(sdof,N);       % displacements vector

% indice vectors pointers
des_id = 1:nidof;
vel_id = nidof+1:2*nidof;
acc_id = 2*nidof+1:3*nidof;
Fob_id = 3*nidof+1:4*nidof;

% separated vectors
des_invertido(iddof,:) = yinv(des_id,:);
vel_invertido(iddof,:) = yinv(vel_id,:);
acc_invertido(iddof,:) = yinv(acc_id,:);
Fob_invertido(iddof,:) = yinv(Fob_id,:);



figure()
% acceleration in mass 1
subplot(2,1,1)
plot(t,acc(3,:),'b')
hold on
plot(t,acc_invertido(3,:),'g--')
xlim([0 50])
ylim([-50 50])
title('Acceleration at massa 1')
xlabel('time')
ylabel('a')

% acceleration in mass 2
subplot(2,1,2)
plot(t,acc(5,:),'b')
hold on
plot(t,acc_invertido(5,:),'g--')
xlim([0 50])
ylim([-50 50])
title('Acceleration at massa 1')
xlabel('time')
ylabel('a')
