%% % dual function theta(z) and its gradient
% Input:
%   G         the given symmetric matrix to be approximated
%   ConstrA:
%        e       the right hand side of equality constraints
%        I_e     row indices of the fixed elements
%        J_e     column indices of the fixed elements
%        l       the right hand side of lower bound constraint
%        I_l     row indices of the lower bound elements
%        J_l     column indices of the lower bound elements
%        u       the right hand side of upper bound constraint
%        I_u     row indices of the upper bound elements
%        J_u     column indices of the upper bound elements
%   n: dimension of the matrix in the primal approximating problem
% Output:
%   theta:
%       dual function value 
%   g_ze,g_zl,g_zu:
%       graidient of dual function correspinding to equlity constrants, lower bounds and upper bounds
%   eig_time, P, lambda:
%        time used for eigen decomposition, matrix composed by eigenvector and the lambda is a vector composed by eigen values

function  [theta,g_ze,g_zl,g_zu,eig_time,P,lambda] = ...
   thetafun(G,e,z_e,I_e,J_e,l,z_l,I_l,J_l,u,z_u,I_u,J_u,n)

k_e = length(e);
k_l = length(l);
k_u = length(u);
k   = k_e + k_l + k_u;

g_ze = zeros(k_e,1);
g_zl = zeros(k_l,1);
g_zu = zeros(k_u,1);

X = zeros(n,n);
for i=1:k_e
  X(I_e(i), J_e(i)) = z_e(i);
end
for i=1:k_l
  X(I_l(i), J_l(i)) = z_l(i) + X(I_l(i), J_l(i));  
end
for i=1:k_u
  X(I_u(i), J_u(i)) = -z_u(i) + X(I_u(i), J_u(i));  %%% upper bound
end
X = 0.5*(X + X');
X = G + X;
X = (X + X')/2;                                     
t1         = clock;
[P,lambda] = MYmexeig(X);
eig_time   = etime(clock,t1);
Xplus=P*diag(max(0,lambda))*P';                      %%%向正部投影
bz=e'*z_e+l'*z_l-u'*z_u;
theta=norm(Xplus,'fro')^2/2-bz;

for i=1:k_e
  g_ze(i)=Xplus(I_e(i), J_e(i));
end
g_ze=g_ze-e;
for i=1:k_l
  g_zl(i)=Xplus(I_l(i), J_l(i));
end
g_zl=g_zl-l;
for i=1:k_u
  g_zu(i)=-Xplus(I_u(i), J_u(i));  
end
g_zu=g_zu+u;
return
