%%  Ref_sub_NewtonStep:
% take the refined adaptive accelerated strategy
% used in AASA(Adaptive_Ref) (AASA_BB_Lbfgs_SNN_adapt_Ref)
%%
function [y_new, g_y_new, NewtonStep_succ, num_Newton, NewtonStep_info, good_cut, theta_y, epsilonk, f_eval, g_eval] = ...
    Ref_sub_intera(G,   e,z_e,I_e,J_e,   l,z_l,I_l,J_l,   u,z_u,I_u,J_u,   n,...
                   zc_lu,BBstepsize,    zb_l,theta,g,P,lambda,   trho,cgtol,cgmaxit,tol,    epsilonk,f_eval,g_eval)
%% detect the active set (C) and free set (F)
em=min(sqrt(epsilonk)*1e-4,1e-6);
% em=min(epsilonk,1e-6);
k_e = length(z_e); k_l = length(z_l); k_u = length(z_u);
k= k_e + k_l + k_u;
z_lu = [z_l;z_u];
%compute the working set C for estimating the active variables at BB point
C=find(abs(zc_lu-zb_l<em) | (g(k_e+1:end)<=-em/2/BBstepsize & abs(z_lu-zb_l<em/2)));
C=C+k_e;
C_old = C;
num_C_old = length(C_old);
fprintf(' #C:%d', num_C_old)
%compute the working set F for estimating the free variables at BB point
F=(1:k)';
F(C)=[];
y = zeros(k,1);
y_temp = [z_e;z_lu];
y(F) = y_temp(F);

% nonsmooth Newton method acceleration
% [theta,g_ze_y,g_zl_y,g_zu_y,eig_time_loc,P,lambda] ...
%      = thetafun(G,e,z_e,I_e,J_e,l,z_l,I_l,J_l,u,z_u,I_u,J_u,n);
[theta_y, g_ye_new, g_yl_new, g_yu_new, eig_time_loc, P, lambda]  ...
         = thetafun(G,e,y(1:k_e),I_e,J_e,...
         l,y(k_e+1:k_e+k_l),I_l,J_l, u,y(end-k_u+1:end),I_u,J_u,n);
    g_ylu_new=[g_yl_new;g_yu_new]; 
    g_y_new = [g_ye_new; g_ylu_new];   % 更新梯度以备下次使用
    g = g_y_new;
    f_eval = f_eval + 1;
    g_eval = g_eval + 1;
    
    
if theta_y <=theta
    good_cut = 1; fprintf('\n =>Good Cut:   0  %.3e  %.3e   #C%d;',theta_y,epsilonk,num_C_old)
else
    good_cut = 0; fprintf('\n =>Bad Cut:    0  %.3e  %.3e   #C%d;',theta_y,epsilonk,num_C_old)
end
theta_y_old = theta_y;
g_ye = g(1:k_e); g_yl = g(k_e+1:k_e+k_l); g_yu = g(end-k_u+1:end);
g_ylu = [g_yl;g_yu];     g_y = [g_ye; g_ylu];


%%
num_Newton = 0;    NewtonStep_succ = 0;     NewtonStep_info = 0;
while 1
    Omega12 = omega_mat2(lambda,n);
    %%%%% 计算搜索方向
    %c = precond_matrix2(k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n); % comment this line for  no preconditioning
    c = ones(k_e+k_l+k_u,1);
    [d,cgflag,cgrelres,iterk]  = pre_cg2(-g_y(F),cgtol,cgmaxit,c(F),F, k_e,k_l,k_u,I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P,n,k);
    ds_y = zeros(k,1);
    ds_y(F) = d;
    ds_e_y = ds_y(1:k_e);   ds_lu_y = ds_y(k_e+1:end);
    
    %%%% 根据搜索方向更新
    y_new = y + ds_y; % 单走牛顿步
    y_new(k_e+1:end) = proje(y_new(k_e+1:end), zb_l); % 走完牛顿步，投影
    y_new_e = y_new(1:k_e);     
    y_new_lu = y_new(k_e+1 : end);
    y_new_l = y_new_lu(1:k_l);  y_new_u = y_new_lu(k_l+1:end);
    deriv = trho*(g_ye'*ds_e_y + g_ylu'*ds_lu_y);
    
    [theta_y, g_ye_new, g_yl_new, g_yu_new, eig_time_loc, P, lambda]  ...
         = thetafun(G,e,y_new_e,I_e,J_e,l,y_new_l,I_l,J_l,u,y_new_u,I_u,J_u,n);
    g_ylu_new=[g_yl_new;g_yu_new]; 
    g_y_new = [g_ye_new; g_ylu_new];   % 更新梯度以备下次使用
    f_eval = f_eval + 1;
    g_eval = g_eval + 1;
    
    num_Newton = num_Newton+1;
    
    % 计算残差
    pg=[-g_ye_new;  proje(y_new_lu-g_ylu_new,zb_l) - y_new_lu];
    epsilonk=norm(pg);
    em=min(sqrt(epsilonk)*1e-4,1e-6);
    C=k_e+find((g_y_new(k_e+1:end)>=0 & abs(y_new_lu-zb_l<em/2)));
    
    
    
    
    % 判断是否停止
    if epsilonk<= tol
        NewtonStep_succ = 1;
        out.msg='converge';
        fprintf('\n Newton step:  %d  %.3e  %.3e  #C:%d; ',num_Newton,theta_y,epsilonk,length(C))
        fprintf('\n total %d Newton Steps are ACCEPTED; Convergens;',num_Newton);
        break;      
    end
    
    if theta_y >= theta_y_old + deriv/100 
        % 如果牛顿步下降量不被接受，则退出牛顿步内循环，返回进入牛顿步循环的点（走LBFGS）
        fprintf('\n Newton step:  %d  %.3e  %.3e  #C:%d; ',num_Newton,theta_y,epsilonk,length(C))
        fprintf('\n Quit due to UNACCEPTED Newton Setp %d is  with nrmg:%.3e  #C:%d; ',num_Newton+1,epsilonk,length(C))
        
        y_new = y;             % 反向更新并输出y_new（上一步的）
        g_y_new = g_y;         % 反向更新并输出梯度（上一步的）
        NewtonStep_info = 1;   % 由于牛顿步不下降而退出
        break;
    else % 如果牛顿步下降量被接受，则接受新点y_new及其梯度g_y_new
        if length(C) ~= num_C_old || sum(setdiff(C, C_old))>0
        fprintf('\n Newton step:  %d  %.3e  %.3e  #C:%d; ',num_Newton,theta_y,epsilonk,length(C))
        fprintf('\n Quit due to Different Active Set. #new_C:#C %d:%d ', length(C),num_C_old)
        NewtonStep_info = 2;    % 由于Active Set不一致而退出
        break;
        end
        fprintf('\n Newton step:  %d  %.3e  %.3e  #C:%d; ',num_Newton,theta_y,epsilonk,length(C))
        theta_y_old = theta_y; % 更新目标函数以备下次使用
        y = y_new;             % 更新y以备下次使用
        g_y = g_y_new;         % 更新梯度以备下次使用
        g_ye = g_ye_new;     g_ylu = g_ylu_new;
    end
end
end



