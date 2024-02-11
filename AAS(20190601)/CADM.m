function [X,Z,info,itr,time_used]=CADM(G,ConstrA,OPTIONS)
%%%%%%%%%%%%% This code is designed to solve %%%%%%%%%%%%%%%%%%%%%
%%       min    0.5*<X-G, X-G>
%%       s.t.   X_ij  = e_ij     for (i,j) in (I_e, J_e)
%%              X_ij >= l_ij     for (i,j) in (I_l, J_l)
%%              X_ij <= u_ij     for (i,j) in (I_u, J_u)
%%              X    >= tau0*I   X is SDP (tau0>=0 and may be zero)
%%
%   Parameters:
%   Input
%   G         the given symmetric matrix
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
%   OPTIONS   parameters in the OPTIONS structure
%   Z0        the initial guess of dual variables
%
%   Output
%   X         the optimal primal solution

%   infos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%-----------------------------------------
%%% get constraints infos from constrA
%%-----------------------------------------
%%
e   = ConstrA.e; I_e = ConstrA.Ie; J_e = ConstrA.Je;
l   = ConstrA.l; I_l = ConstrA.Il; J_l = ConstrA.Jl;
u   = ConstrA.u; I_u = ConstrA.Iu; J_u = ConstrA.Ju;
k_e = length(e); k_l = length(l);  k_u = length(u);
k   = k_e + k_l + k_u;  n = size(G,1);

%%
%%-----------------------------------------
%% get parameters from the OPTIONS structure.
%%-----------------------------------------
%%

if exist('OPTIONS')
    if isfield(OPTIONS,'tol');         tol        = OPTIONS.tol;    else,    tol     = 1.0e-6; end
    if isfield(OPTIONS,'maxit');       maxit      = OPTIONS.maxit;  else,    maxit   = 3000;   end
    if isfield(OPTIONS,'disp');        DISP       = OPTIONS.disp;   else,    DISP    = 1;      end
    if isfield(OPTIONS,'maxtime');     maxtime    = OPTIONS.maxtime;else,    maxtime = 1500;   end
end

%%%%----------------------------------------------------------------------
if DISP
    fprintf('\n ******************************************************** \n')
    fprintf( '    The CADM method       ')
    fprintf('\n ******************************************************** \n')
    fprintf('\n The information of this problem is as follows: \n')
    fprintf(' Dim. of    sdp      constr  = %d \n',length(G))
    fprintf(' Num. of equality    constr  = %d \n',k_e)
    fprintf(' Num. of lower bound constr  = %d \n',k_l)
    fprintf(' Num. of upper bound constr  = %d \n',k_u)
    fprintf(' The lower bounds: [ %2.1f, %2.1f ] \n',min(l),max(l))
    fprintf(' The upper bounds: [ %2.1f, %2.1f ] \n',min(u),max(u))
end
%%%%----------------------------------------------------------------------


t0 = clock;
k=0;
beta=1;
gammaY=1.5;
Y=eye(n);
Z=zeros(n);
diffYZ=inf;
while diffYZ>tol && k<maxit
    Y0=Y; Z0=Z; k=k+1;
    X=(Y0*beta+Z0+G)/(1+beta);
    [V,D]=mexeig(X); D=max(0,D); X=V*diag(D)*V';
    Y=(X*beta-Z0+G)/(1+beta);
    for i=1:k_e
        Y(I_e(i), J_e(i)) = e(i);
    end
    for i=1:k_l
        Y(I_l(i), J_l(i)) = max( Y(I_l(i), J_l(i)),l(i));
    end
    for i=1:k_u
        Y(I_u(i), J_u(i)) = min( Y(I_u(i), J_u(i)),u(i));
    end
    Y=Y-tril(Y,-1)+triu(Y,1)';
    %%%%%%%%% Calculating the step size %%%%%%%%%%%%%%
    EY=Y0-Y; EZ=(X-Y)*beta;
    T1 = EY(:)'*EY(:); T2 = EZ(:)'*EZ(:); TA=T1*beta + T2/beta;
    T2 = (EY(:)'*EZ(:));
    alpha=(TA-T2)*gammaY/TA;
    Y=Y0-EY*alpha; Z=Z0-EZ*alpha;
    
    diffYZ=max(max(max(abs(Y-Y0))),max(max(abs(Z-Z0))));
    tt = etime(clock,t0);
    [hh,mm,ss] = time(tt);
    fprintf('k=%4d epsm=%9.3e %d:%d:%d \n',k,diffYZ,hh,mm,ss);
    % verify max time_used
    if tt>maxtime
        fprintf('\n *********************************************************************');
        fprintf('\n *********************************************************************');
        fprintf('\n Warning: Progress is too slow! Stop if the follows happend :');
        fprintf('\n Total computing time is beyond Max time: %ds',maxtime);
        fprintf('\n *********************************************************************')
        fprintf('\n *********************************************************************');
        break;
    end
    
end
itr=k;
if itr<=maxit
    info=-1;
else
    info=1;
end
time_used = etime(clock,t0);
%%% To change the format of time
function [h,m,s] = time(t)
t = round(t);
h = floor(t/3600);
m = floor(rem(t,3600)/60);
s = rem(rem(t,60),60);
%%% End of time.m
