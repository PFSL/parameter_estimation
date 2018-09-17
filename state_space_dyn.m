function [y,x] = state_space_dyn(nidof, N, dt, Mef, Kef, Cef, Fef)
%% solution at state space

%% mount state matrix
Ac = [
    [zeros(nidof,nidof), eye(nidof,nidof)];
    [-Mef\Kef, -Mef\Cef]
    ];

%% mount input matrix
Bc = [
    [zeros(nidof,nidof)];
    [inv(Mef)]
    ];

%% mount observation matrix
% order to observe is: displacement | velocities | accelerations | forces
C = [
    [eye(nidof,nidof), zeros(nidof,nidof)];
    [zeros(nidof,nidof), eye(nidof,nidof)];
    [-Mef\Kef, -Mef\Cef];
    [zeros(nidof,nidof), zeros(nidof,nidof)]
    ];

%% mount direct transmission matrix
% order to transmission is: displacement | velocities | accelerations | forces
D = [
    [zeros(nidof,nidof)];
    [zeros(nidof,nidof)];
    [inv(Mef)];
    [eye(nidof,nidof)]
    ];

% integration procedure
A = expm(Ac*dt);
B = Ac\(A-eye(2*nidof))*Bc;

% taken effective forces applyed
u = Fef;

% define initial state (bc equals to zero in displacements and velocities)
x = zeros(2*nidof,1);

%% time march procedure (linear acceleration)
for k = 1:N
    y(:,k) = C*x + D*u(:,k);
    x = A*x + B*u(:,k);
end
