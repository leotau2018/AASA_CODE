function [X,z,info,Iterations,time_used,condn] = AASA_BB_Lbfgs_SNN_adapt_Ref(G,ConstrA,OPTIONS,z0)
%AASA_BB_Lbfgs_SNN_adapt_Ref，测试的主程序。主要思路：先利用了BB-point的信息，将变量分成积极约束集和非积极约束集，然后返回x_k处迭代。
%积极约束集直接拉向边界，非积极约束集采用LBFGS法或半光滑牛顿法迭代。
%同时考虑了长短BBstep size的选择。
%Refined adaptive acceleration strategy for Algotithm 1
%%%%%%%%%%%%% This code is designed to solve %%%%%%%%%%%%%%%%%%%%%
%%       min    0.5*sigma_1*<X-G, X-G> + sigma_2<X>
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
%   z0        the initial guess of dual variables
%
%   Output
%   X         the optimal primal solution
%   z:
%      z.e    the optimal dual solution to equality constraints
%      z.l    the optimal dual solution to lower bound constraints
%      z.u    the optimal dual solution to upper bound constraints
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
k   = k_e + k_l + k_u;  n = length(G);

%%
%%-----------------------------------------
%% get parameters from the OPTIONS structure.
%%-----------------------------------------
%%

if exist('OPTIONS','var')
    if isfield(OPTIONS,'tol');             tol      = OPTIONS.tol;     else,  tol     = 1.0e-6;  end
    if isfield(OPTIONS,'maxit');           maxit    = OPTIONS.maxit;   else,  maxit   = 2000;  end
    if isfield(OPTIONS,'cgtol');           cgtol    = OPTIONS.cgtol;   else,  cgtol   = 1.0e-2; end
    if isfield(OPTIONS,'cgmaxit');         cgmaxit  = OPTIONS.cgmaxit; else,  cgmaxit = 200;   end
    if isfield(OPTIONS,'M');               M        = OPTIONS.M;       else,  M       =  5;    end
    if isfield(OPTIONS,'disp');            disp     = OPTIONS.disp;    else,  disp    = 1; end
    if isfield(OPTIONS,'tau');             tau    = OPTIONS.tau; else, tau       = 1;   end
    if isfield(OPTIONS,'eta');             eta    = OPTIONS.eta; else, eta       = 0.2; end
    if isfield(OPTIONS,'Cumpt_Condn');     Cumpt_Condn  = OPTIONS.Cumpt_Condn; else, Cumpt_Condn = 0; end
    if isfield(OPTIONS,'maxtime');         maxtime      = OPTIONS.maxtime; else, maxtime=1500;end
    if isfield(OPTIONS,'Ada_on');          Ada_on = OPTIONS.Ada_on; else, Ada_on = 1;   end
    if isfield(OPTIONS,'Ref_on');          Ref_on = OPTIONS.Ref_on; else, Ref_on = 1;   end
    if isfield(OPTIONS,'e_1');             e_1    = OPTIONS.e_2;   else, e_1     = 2e-1;   end
    if isfield(OPTIONS,'e_2');             e_2    = OPTIONS.e_2;   else, e_2     = 5e-3;   end
end


t0 = clock;

%%% reset the pars
for i = 1:k_e
    G(I_e(i),J_e(i)) = e(i);
    if I_e(i) ~= J_e(i)
        G(J_e(i),I_e(i)) = e(i);
    end
end
%%%%%加入不等式约束条件

%%%
f_eval    = 0;
g_eval    = 0;
eig_time  = 0;
theta_time  = 0;


start_index = 1;
S = zeros(k,M);
Y = S;
RHO=zeros(M,1);
lr=1;
gamma = 1;gammaL=1;
trho=1e-2;
epsilonk=inf;
hist_epsilon=[]; %history of KKT error
hist_rhol=[inf]; % history of contribution ratios of L-BFGS
hist_rhon=[inf]; % history of contribution ratios of semi-smooth Newton
cons_lbfgs=0;
cons_ssn=0; cons_ssn_temp = 0; NewtonStep_Info = 0; num_Newton=0;
ratios=1;
ratios1=1;
acc_type=1;% acceleration type
%%% initial value

if ( nargin == 4 )
    z_e = z0.e;
    z_l = z0.l;
    z_u = z0.u;
else
    z_e = zeros(k_e,1);
    %z_l = zeros(k_l,1);
    %z_u = zeros(k_u,1);
    z_l=0.1*ones(k_l,1);
    z_u=0.1*ones(k_u,1);
end

%%%%对偶问题初始点赋值。
z0_e = z_e;
z0_l = z_l;
z0_u = z_u;
z0_lu =[z0_l;z0_u];

if disp
    fprintf('\n ******************************************************** \n')
    if Ada_on == 1 
        if Ref_on ==1, fprintf( '    The Redined AASA(adaptive) method       ')
        else,          fprintf( '    The AASA(adaptive) method       ')
        end
    else,              fprintf( '    The AASA(L-BFGS) method       ')
    end
    fprintf('\n ******************************************************** \n')
    fprintf('\n The information of this problem is as follows: \n')
    fprintf(' Dim. of    sdp      constr  = %d \n',n)
    fprintf(' Num. of equality    constr  = %d \n',k_e)
    fprintf(' Num. of lower bound constr  = %d \n',k_l)
    fprintf(' Num. of upper bound constr  = %d \n',k_u)
    fprintf(' The lower bounds: [ %2.1f, %2.1f ] \n',min(l),max(l))
    fprintf(' The upper bounds: [ %2.1f, %2.1f ] \n',min(u),max(u))
end


t1 = clock;
[theta,g_ze,g_zl,g_zu,eig_time_loc,P,lambda] = thetafun(G,e,z_e,I_e,J_e,l,z_l,I_l,J_l,u,z_u,I_u,J_u,n);
theta_time = theta_time + etime(clock,t1);
eig_time   = eig_time + eig_time_loc;
f_eval = f_eval + 1;
g_eval = g_eval + 1;
Q = 1; Cval = theta; tgamma = .85;
zb_l=zeros(k_l+k_u,1);                %%%原问题的等式约束对应的对偶变量无约束，原问题的不等式约束对应的对偶变量只有下约束，记为zb_l
itr=1;
if disp
    
    tt = etime(clock,t0);
    [hh,mm,ss] = time(tt);
    fprintf('\n    Iter    theta       nrmg         StepLen  time_used type #C ')
    fprintf('\n    %2.0f   %3.2e     %3.2e      %3.2e  %d:%d:%d ',0,theta,epsilonk,tau,hh,mm,ss)
end

%% main iteration
for itr = 1 : maxit
    gp_ze=g_ze;
    gp_zl=g_zl;
    gp_zu=g_zu;
    gp_zlu=[gp_zl;gp_zu];
    g_zlu =[g_zl;g_zu];
    k_lu  =k_l+k_u;       %k_lu是不等式约束对应的乘子个数，对应的，k_e是等式约束对应的乘子个数
    
    %%%%%%% Compute the BB point %%%%%%%%%%%%%%%%%%%%%%%%
    
    if (gamma/gammaL<0.2)   %选择自适应步长
        BBstepsize=gamma;
    else
        BBstepsize=gammaL;
    end
    
    BBstepsize = max(min(BBstepsize, 1e20), 1e-20);
    
    zc_lu=proje(z0_lu-BBstepsize*g_zlu,zb_l); %zc_lu 是BB point
    g=[g_ze;g_zl;g_zu]; %综合成一个向量。
    
    if Ada_on == 1
        if epsilonk<=e_1*ratios % e_1 = 2e-1
            acc_type=2;
            cons_ssn_temp = cons_ssn_temp+1;
            if epsilonk>=5e-2
                cgmaxit=10;
            elseif epsilonk>=2e-2
                cgmaxit=15;
            elseif epsilonk>=1e-3
                cgmaxit=20;
            elseif epsilonk>=1e-4
                cgmaxit=25;
            end
        else
            acc_type=1; bad_Newton = 0;
        end

        if acc_type==2 && (tau<1 || chi<.3) &&  rem(cons_ssn_temp,3)==0
            if epsilonk>=5e-2 || bad_Newton >= 1
                acc_type=1;
                cons_ssn_temp=0;
                cons_lbfgs=cons_lbfgs+1;
                bad_Newton = 0;
            else
                bad_Newton = bad_Newton+1;
            end
        end
        %% take the refined adaptive accelerated strategy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if acc_type == 2 && epsilonk < e_2*ratios1 && Ref_on  % e_2 = 5e-3
            cal_Newtonstep_time0 = clock;
            [y_new, g_y_new, NewtonStep_succ, num_Newton, NewtonStep_Info(itr), good_cut, theta_y, epsilonk, f_eval, g_eval] = ...
                Ref_sub_intera(G,   e,z_e,I_e,J_e,   l,z_l,I_l,J_l,   u,z_u,I_u,J_u,   n,...
                zc_lu,BBstepsize,    zb_l,theta,g,P,lambda,   trho,cgtol,2000,tol,    epsilonk,f_eval,g_eval);
            cal_Newtonstep_time = etime(clock,cal_Newtonstep_time0);
            fprintf(' Time:%.2fs', cal_Newtonstep_time)
            cons_ssn=cons_ssn+num_Newton;
            if NewtonStep_succ == 1
                z_e = y_new(1:k_e);  % 计算成功，输出变量，用于后续计算矩阵X
                z_l = y_new(k_e+1:k_e+k_l);
                z_u = y_new(end-k_u+1:end);
                %z_lu = y_new(k_e+1:end);
                %g = g_y_new;
                if disp
                    tt = etime(clock,t0);
                    [hh,mm,ss] = time(tt);
                    fprintf('\n == %2.0f   %3.5e   %3.5e   %3.2e  %d:%d:%d  tp:Newton Step',itr,theta_y,epsilonk,tau,hh,mm,ss);
                end
                out.msg='converge';
                fprintf('\n It converges!!!');
                break;
            else
                if num_Newton >= 1 || good_cut
                    % 内层牛顿循环未达到最优点，变量更新为内层牛顿循环得到的可行点
                    % 更新点为下降点，好的cuted点或上一次下降的牛顿步点
                    z0_e = y_new(1:k_e); 
                    % z0_l = y_new(k_e+1:k_e+k_l);
                    % z0_u = y_new(end-k_u+1:end);
                    z0_lu = y_new(k_e+1:end);
                    g = g_y_new;
                    g_ze = g_y_new(1:k_e);  g_zlu = g_y_new(k_e+1:end);
                    acc_type = 1; % or acc_type = 2
                    ratios1=ratios1/5;
                else
                    % 内层牛顿循环未达到最优点，且 cut不好并且第一次牛顿步不下降
                    acc_type = 1; % or acc_type = 2
                    ratios1=ratios1/5;
                end
            end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if acc_type == 2
            cons_ssn=cons_ssn+1;
        end
    else, num_Newton = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%% get the search direction %%%%%%%%%%%%%%%%
    calstep_time0 = clock;
    [ds_e, ds_lu, C] = getD_byBB_new(S, Y, lr-1, start_index,  gamma, M, P, lambda,n,k_e,k_l,k_u, I_e,J_e, I_l,J_l,I_u, J_u, zc_lu, zb_l,epsilonk,BBstepsize, z0_lu,g,cgtol,cgmaxit,acc_type);%get the search direction at current poin，w是omega.
    %[ds_e,ds_lu]=getD_byCauchy( S, Y, lr-1, start_index, gamma, M,k_e,k_lu, zc_lu, zb_l,epsilonk, z0_lu,g);
    calstep_time = etime(clock,calstep_time0);
    %%%%search direction completed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tau=1;
    nls=1;
    while 1
        % calculate g, theta,
        
        z_e =  z0_e + tau*ds_e;                %tau=1
        z_lu = z0_lu + tau*ds_lu;              %tau=1,此时走一步已经使得部分（C2中）分量走到l-bond上了
        z_lu = proje(z_lu,zb_l);
        z_l = z_lu(1:k_l);     z_u=z_lu(k_l+1:end);
        deriv = trho*(g_ze'*ds_e+g_zlu'*ds_lu);
        
        
        t1 = clock;
        [theta,g_ze,g_zl,g_zu,eig_time_loc,P,lambda] = thetafun(G,e,z_e,I_e,J_e,l,z_l,I_l,J_l,u,z_u,I_u,J_u,n);
        g_zlu=[g_zl;g_zu];
        theta_time = theta_time + etime(clock,t1);
        eig_time   = eig_time + eig_time_loc;
        f_eval = f_eval + 1;
        g_eval = g_eval + 1;
        
        if theta <= Cval + tau*deriv || nls >= 10
            break;
        end
        tau = eta*tau;
        if acc_type==2 && tau<1e-3
            acc_type=1;
            tau=1;
            nls=1;
            ratios=ratios/10;
            [ds_e, ds_lu] = getD_byBB_new(S, Y, lr-1, start_index,  gamma, M, P, lambda,n,k_e,k_l,k_u, I_e,J_e, I_l,J_l,I_u, J_u, zc_lu, zb_l,epsilonk,BBstepsize, z0_lu,g,cgtol,cgmaxit,acc_type);
            continue;
        end
        nls = nls+1;
    end
    %%%终止条件
    %pg=[dNLS_e;dNLS_lu];
    pg=[-g_ze;proje(z_lu-g_zlu,zb_l)-z_lu];
    epsilonk=norm(pg);
    hist_epsilon=[hist_epsilon;epsilonk];
    if itr>1 && cons_ssn>=4
        chi=(hist_epsilon(end-3)-epsilonk)/epsilonk;
    else
        chi=inf;
    end
    %acc_type
    if epsilonk<= tol
        
        if disp
            tt = etime(clock,t0);
            [hh,mm,ss] = time(tt);
            fprintf('\n    %2.0f   %3.5e   %3.5e   %3.2e  %d:%d:%d  tp:%d',itr,theta,epsilonk,tau,hh,mm,ss,acc_type);
        end
        out.msg='converge';
        fprintf('\n It converges!!!');
        break;
    end
    s =  [z_e;z_lu]- [z0_e;z0_lu];
    y =  [g_ze-gp_ze; g_zlu-gp_zlu];
    sy = s'*y;
    if sy<0
        sy=-sy;
        y=-y;
    end
    ss = s'*s;
    yy = y'*y;
    Qp = Q;    Q = tgamma*Qp + 1; Cval = (tgamma*Qp*Cval + theta)/Q;
    
    
    if(ss>1e-12 && sy>eps)
        rho = 1 / sy;
        gamma = sy / yy;%这里的gamma其实相当于短BB-step size
        %计算long BB-step size
        gammaL=ss/sy;   %记gammaL为为长BB-step size
        
        
        %% transport discard and store new.
        if(lr <= M)
            Y(:,lr) = y;
            S(:,lr)= s;
            RHO(lr) = rho;
            lr = lr+1;
        else
            Y(:,start_index) = y;
            S(:,start_index) = s;
            RHO(start_index) = rho;
            start_index = start_index + 1;
            if(start_index > M)
                start_index = 1;
            end
        end
        
    end
    z0_e = z_e;
    z0_l = z_l;
    z0_u = z_u;
    z0_lu = z_lu;
    
    if disp
        tt = etime(clock,t0);
        [hh,mm,ss] = time(tt);
        fprintf('\n    %2.0f   %3.5e   %3.5e   %3.2e  %d:%d:%d  tp:%d',itr,theta,epsilonk,tau,hh,mm,ss,acc_type);
    end
    if tt>maxtime
        fprintf('\n *********************************************************************');
        fprintf('\n *********************************************************************');
        fprintf('\n Warning: Progress is too slow! Stop if the follows happend :');
        fprintf('\n Total computing time is beyond Max time: %ds',maxtime);
        fprintf('\n *********************************************************************')
        fprintf('\n *********************************************************************');
        break;
    end
    if itr>=maxit    %异常情况下记录该案例
        mkdir('Log')
        save Log/oddOPTIONS.mat OPTIONS ;
        save Log/oddConstrA.mat ConstrA ;
        save Log/oddG.mat       G ;
        fprintf('\n *********************************************************************');
        fprintf('\n Warning: Progress is too slow! Stop if the follows happend :');
        fprintf('\n Total computing time is beyond Max iteration: %d',maxit);
        fprintf('\n *********************************************************************')
    end
end

% optimal primal solution X*
X = zeros(n,n);
for i=1:k_e
    X(I_e(i), J_e(i)) = z_e(i);                       %%% equality constraint
end
for i=1:k_l
    X(I_l(i), J_l(i)) = z_l(i) + X(I_l(i), J_l(i));   %%% lower bound
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
X=P*diag(max(0,lambda))*P';

z.e = z_e;
z.l = z_l;
z.u = z_u;
info=1;
if disp==1
    time_used = etime(clock,t0);
    %fid = fopen('result.txt','wt');
    %fprintf(fid,'\n');
    fprintf('\n\n ================ Final Information ================= \n');
    fprintf(' Total number of iterations      = %2.0f \n', itr);
    fprintf(' Number of func. evaluations     = %2.0f \n', f_eval);
    fprintf(' Number of grad. evaluations     = %2.0f \n', g_eval);
    %     fprintf(' Primal objective value          = %d \n', prim_val);
    %     fprintf(' Dual objective value            = %d \n', -dual_val);
    fprintf(' Computing time for eigen-decom        = %3.1f \n', eig_time);
    fprintf(' Computing time for the merit fun.     = %3.1f \n', theta_time);
    fprintf(' Total computing time (secs)           = %3.1f \n', time_used);
    fprintf(' ====================================================== \n');
    %fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Cumpt_Condn == 1
fprintf('\n\n ================ Conditon Number Computing ================= \n');
cond_time0=clock;
condn=Condnum(k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u,P, n, lambda, zc_lu, zb_l,epsilonk,BBstepsize, z0_lu,g);
condest_time=etime(clock,cond_time0);
fprintf(' Conditon Number of the optimial dual Hessian   = %.4e \n', condn);
fprintf(' Extra time for computing Condition Number      = %3.1fs \n', condest_time);
fprintf(' ==========================END=============================== \n');
else
condn = 0;
end
Iterations.total_itr = itr-1+num_Newton;
Iterations.cons_ssn = cons_ssn;
Iterations.NewtonStep_BreakbyUnacceptedReduce = sum(NewtonStep_Info==1); % 由于牛顿步不下降而退出
Iterations.NewtonStep_BreakbyDiffActiveSet = sum(NewtonStep_Info==2);    % 由于Active Set不一致而退出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% end of the main program %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  **************************************
%%  ******** All Sub-routines  ***********
%%  **************************************

%%% To change the format of time
function [h,m,s] = time(t)
t = round(t);
h = floor(t/3600);
m = floor(rem(t,3600)/60);
s = rem(rem(t,60),60);
%%% End of time.m






