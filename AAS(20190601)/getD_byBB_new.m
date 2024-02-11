%%compute the search direction based on the information at BB point
% used in AASA(Adaptive) (AASA_BB_Lbfgs_SNN_adapt)
function [ds_e,ds_lu, C]=getD_byBB_new( S, Y, lenS, start_index, gamma, M, P, lambda,n,k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u, zc_lu, zb_l,epsilonk,BBstepsize, z0_lu,g,cgtol,cgmaxit,type)
%type=1 L-BFGS acceleration
%     2 nonsmooth Newton method acceleration
em=min(sqrt(epsilonk)*1e-4,1e-6);
% em=min(epsilonk,1e-6);
k=size(S,1);
ds=[zeros(k_e,1);zc_lu-z0_lu];
%compute the working set C for estimating the active variables at BB point
C=find(abs(zc_lu-zb_l)<em | (g(k_e+1:end)<=-em/2/BBstepsize & abs(z0_lu-zb_l)<em/2));
% C=find(g(k_e+1:end)>=0 & abs(z0_lu-zb_l<em/2) );
% C=find(abs(zc_lu-zb_l<em));
C=C+k_e;
fprintf(' #C:%d', length(C))
%compute the working set F for estimating the free variables at BB point
F=(1:k)';
F(C)=[];

if type==2
    % nonsmooth Newton method acceleration
    Omega12 = omega_mat2(lambda,n);
    %c = precond_matrix2(k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n); % comment this line for  no preconditioning
    c = ones(k_e+k_l+k_u,1);    
    [d,cgflag,cgrelres,iterk]  = pre_cg2(-g(F),cgtol,cgmaxit,c(F),F, k_e,k_l,k_u,I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P,n,k);
    %iterk
    ds(F)=d;
    if cgflag<0
        type=1;
    end
end
if type==1
    % L-BFGS acceleration
    if lenS==M && start_index>1
        ind=[start_index:M, 1:start_index-1];
        S=S(:,ind);
        Y=Y(:,ind);
    else
        ind=1:lenS;
        S=S(:,ind);
        Y=Y(:,ind);
    end
    
    Wkb=[gamma*Y,S];                                                      %LBFGS相关变量计算
    SY=S'*Y;
    Dk=diag(diag(SY));
    Rk=triu(SY);Rk=Rk\eye(lenS);
    Mkb=[zeros(lenS),-Rk;-Rk',Rk'*(Dk+gamma*(Y'*Y))*Rk];
    NNg=zeros(k,1);
    NNg(F)=-g(F);
    tmp=Wkb*(Mkb*(Wkb'*NNg));
    tmp=tmp+gamma*NNg;
    ds(F)=tmp(F);        %ds是假设没有约束时使用拟牛顿法的搜索方向
end

%以下为了格式统一而做调整。
ds_e=ds(1:k_e);           %length of ds_e=k_e
ds_lu=ds(k_e+1:end);        %length of ds_lu=k_lu,此处为约束乘子的非积极约束部分采用拟牛顿方向
end