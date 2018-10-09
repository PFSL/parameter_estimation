function [k_loc, m_loc] = stiff_local2d(L, A, E, rho, beta, ipt)

%% mount local stiffness and mass matrix based on FEM

% director sines and cosines
c = cos(beta);
s = sin(beta);

% sitffness matrix for truss element
k_loc = (A*E/L)*[c*c   c*s   -c*c  -c*s;
             c*s   s*s   -c*s  -s*s;
             -c*c  -c*s   c*c   c*s;
             -c*s  -s*s   c*s   s*s];

% massa matrix for truss element
if ipt == 1
    m_loc = (rho*A*L/6)*[2  0  1  0;    % consistent
                     0  2  0  1;
                     1  0  2  0;
                     0  1  0  2];
                 
else
    m_loc = (rho*A*L/2)*eye(4);         % lumped
end
end