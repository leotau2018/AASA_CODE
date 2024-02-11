
[G,ConstrA]=Case_test(2,2000,[50,50,50]*10);
%% Real Case 
[G,ConstrA]=case_real_test(1,3,792);          % Real Case E7
[G,ConstrA]=case_real_test(2,2,[400,400,400]); % Real Case E8
n=2000;nr=500;ne=n;
total=2*((n-1-nr)*nr+0.5*(nr+1)*nr) +ne
% tau = .0e-8;
% OPTIONS.tau = tau;
OPTIONS.tol = 1e-8;
[X,Z,info,itr,time_used]=CADM(G,ConstrA,OPTIONS);  %CADA
[X0,z0,info,newton_itr,f_eval,newton_time_used] = CaliMat1Mex(G,ConstrA,OPTIONS);  %SNN 半光滑牛顿法
% [X1,z1,F,lbfgs_itr,lbfgs_time_used]      = CaliMatLbfgs2D(G,ConstrA,OPTIONS);  %两搜索方向（樊）
[X,z,info,BB_lbfgs_itr,BB_lbfgs_time_used] =  AASA_BB_Lbfgs_SNN_adapt1(G,ConstrA,OPTIONS); %AASA（l-bfgs）
[X1,z1,F, BB_lbfgsssn_itr,BB_lbfgsssn_time_used]= AASA_BB_Lbfgs_SNN_adapt(G,ConstrA,OPTIONS); %AASA（Adaptive）
% Ref_on = 1(default), means using Refined alforithm based on Alg. 1.
OPTIONS.Ref_on = 1;
[X2,z2,F, BB_lbfgsssn_ref_itr,BB_lbfgsssn_ref_time_used]= AASA_BB_Lbfgs_SNN_adapt_Ref(G,ConstrA,OPTIONS); %AASA（Refined Adaptive）
% ada_on = 0(default 1), means no adaptive strategy used in computing the direction (only Lbfgs for free set)
OPTIONS.Ada_on = 0 ; 
[X2,z2,F, BB_lbfgsssn_ref_itr,BB_lbfgsssn_ref_time_used]= AASA_BB_Lbfgs_SNN_adapt_Ref(G,ConstrA,OPTIONS); %AASA（Refined Adaptive）
Omega12 = omega_mat2(lambda,n);
condn=Condnum(k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n);

[pbfgs_Z,pbfgs_itr, pbfgs_histout,costdata,pbfgs_time_used, pbfgs_flag] = projbfgs(G,ConstrA,OPTIONS);


%% 查看hessian矩阵  normest1()
k=100;cd=3;
V=diag(cd*rand(k,1)+1);
U=orth(rand(k,k));
A=U*V*U^(-1);