%% real example
% 
% [G,ConstrA]=case_real_test(1,3,300);
% [G,ConstrA]=case_real_test(2,2,[300,300,300]);
% In our paper, we tested two cases：
%   (1): [G,ConstrA]=case_real_test(1,3,792);
%        Shanghai Stock Exchange(2016.9-2018.9), n =792, 
%        only lower bounds on all the entries: ll=n，l = -0.1r+0.6
%        denoted by E7 in the paper.
%   (2): [G,ConstrA]=case_real_test(2,2,[400,400,400);
%        Shenzhen Stocks Exchange(2016.9-2018.9), n=1187,
%        400 random positions for lower and upper bounds: e=1,l=--0.3r+1,u=0.3r+0.1。
%        denoted by E8 in the paper.
%
function [G,ConstrA]=case_real_test(realtypeG,casetype,constr_n)
%% get the path
% 设置绝对路径
% datapath = 'D:\wangyunlong\基本材料\studysource\研究生project\3.5 AASA（Accelerated Adaptive Active Set Algorithm）';
% 设置相对路径
datapath = '..';
%% Generate G
if realtypeG==1
%%%%%%%%%%%%%%%%%%%%%% n=792
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
data1=xlsread([datapath,'\Real Examples\cormat.csv']);
G=data1;
n=size(G,1);
%% constr_n为每行约束位置个数。
%当casetype=1时，constr_n为次对角线条数（带状对角宽度），e=1,l=-0.025r+0.5,u=+0.025r+0.5。
%当casetype=2时，constr_n为3*1向量，分别对应le,ll,lu; e=1,l=-0.025r+0.7,u=0.025r+0.7。
%当casetype=3时，只有下约束，ll=n，下约束为-0.1r+0.6

%% Case I  %带状上对角矩阵
%等式约束只在对角线上，且为1,（这个约束在这里无效）；
%不等式约束在次对角线上。上约束和下约束个数相等，等于次对角线条数（带宽）
%下约束为-0.025r+0.5，上约束为+0.025r+0.5
if casetype==1
fprintf('\n Real Test Example,TYPE= RE1 \n')
bread_width=constr_n(1);
%% I_e,J_e
%%%% for fixed  diagonal entries
I_e=[1:1:n]';
J_e=I_e;
k_e=length(I_e);
%% I_l,J_l
%%%  entries with lower bounds
I_l = [];
J_l = [];
for i = 1:n-bread_width
    I_l=[I_l;i*ones(bread_width,1)];
    J_l=[J_l;(i+1:i+bread_width)'];
end
for i = ((n-bread_width)+1):(n-1)
    I_l = [I_l; i*ones(n-i,1)];
    J_l = [J_l;[(i+1):n]'];
end
k_l = length(I_l);

%% I_u,J_u
%%%  entries with lower bounds
I_u = [];
J_u = [];
for i = 1:n-bread_width
    I_u=[I_u;i*ones(bread_width,1)];
    J_u=[J_u;(i+1:i+bread_width)'];
end
for i = ((n-bread_width)+1):(n-1)
    I_u = [I_u; i*ones(n-i,1)];
    J_u = [J_u;((i+1):n)'];
end
k_u = length(I_u);
%% to generate the bound e,l & u
%%%%%%% e
rhs    = ones(n,1);  % diagonal elements
alpha0 = 1;         %通过调整alpha0的值来选择例子中（a）（b）两种情况，alpha0==1时对应情况为（a）
rhs    = alpha0*rhs + (1-alpha0)*rand(n,1);
e      = rhs;
%%%%%%% l
l = -0.025*rand(k_l,1)+0.5;
%%%%%%% u
u = 0.025*rand(k_u,1)+0.5;

max_l = max(l);
min_l = min(l);
max_u = max(u);
min_u = min(u);
end
%% Case II  %随机位置约束
%约束位置随机
%下约束为-0.025r+0.7，上约束为0.025r+0.7
if casetype==2
    fprintf('\n Real Test Example,TYPE= RE2 \n')
    if length(constr_n)<3, constr_n = repmat(constr_n(1),3,1);end
    lh =  constr_n(1);    % number of fixed off diagonal elements in each row
    ll =  constr_n(2);    % number of off diagonal elements of lower bounds in each row
    lu =  constr_n(3);    % number of off diagonal elements of upper bounds in each row
       
    ll = min(ll,n-1);
    lu = min(lu,n-1);
    lh = min(lh,n-1);
    
%     E = 2.0*rand(n,n) - ones(n,n);
%     E = triu(E) + triu(E,1)';
%     
%     E=E'*E;
%     E=E-diag(diag(E))+eye(n);
%     G = E;
    %% I_e,J_e
    %%%% for fixed  diagonal entries
    I_e = [1:1:n]';
    J_e = I_e;
    k_e = length(I_e);
    
    %%  I_l,J_l
    %%%  entries with lower bounds
    I_l = [];
    J_l = [];
    for i = 1:n-ll
        r = rand(n-i,1);
        [r,ind] = sort(r);
        I_l = [I_l; i*ones(ll,1)];
        J_l = [J_l; i+ind(n-i-ll+1:n-i)];
    end
    for i = ((n-ll)+1):(n-1)
        I_l = [I_l; i*ones(n-i,1)];
        J_l = [J_l;[(i+1):n]'];
    end
    k_l = length(I_l);
    %%  I_u,J_u
    %%%%%  entries with upper bounds
    I_u = [];
    J_u = [];
    for i = 1:n-lu
        r = rand(n-i,1);
        [r,ind] = sort(r);
        I_u = [I_u; i*ones(lu,1)];
        J_u = [J_u; i+ind(n-i-lu+1:n-i)];
    end
    for i = ((n-lu)+1):(n-1)
        I_u = [I_u; i*ones(n-i,1)];
        J_u = [J_u;[(i+1):n]'];
    end
    k_u = length(I_u) ;
    %% to generate the bound e,l & u
    %%%%%%% e
    rhs    = ones(n,1);  % diagonal elements
    alpha0 = 1;          %通过调整alpha0的值来选择例子中（a）（b）两种情况，alpha0==1时对应情况为（a）
    rhs    = alpha0*rhs + (1-alpha0)*rand(n,1);
    e      = rhs;
    %%%%%%% l
    l = -0.025*rand( k_l,1)+0.7;
    %l = 0.50 * (2*rand(k_l,1)-ones(k_l,1));
    %l = 1.0 * (rand(k_l,1) - ones(k_l,1));
    
    %%%%%%% u
    u = 0.025*rand(k_u,1)+0.7;
    %u = 1.0*(rand(k_l,1) - ones(k_l,1));
    max_l = max(l);
    min_l = min(l);
    max_u = max(u);
    min_u = min(u);
    
end
%% Case III  %只有下约束且ll=n
%下约束为-0.1r+0.6
if casetype==3
    fprintf('\n Real Test Example,TYPE=  RE3 \n')
    fprintf('\n 只有下约束-0.1r+0.6,所有位置%d \n', n)
    
%     lh =  constr_n(1);    % number of fixed off diagonal elements in each row
%     ll =  constr_n(2);    % number of off diagonal elements of lower bounds in each row
%     lu =  constr_n(3);    % number of off diagonal elements of upper bounds in each row
     
        
    lh =  0;    % number of fixed off diagonal elements in each row
    ll =  n;    % number of off diagonal elements of lower bounds in each row
    lu =  0;    % number of off diagonal elements of upper bounds in each row
    
    ll = min(ll,n-1);
    lu = min(lu,n-1);
    lh = min(lh,n-1);
    
%     E = 2.0*rand(n,n) - ones(n,n);
%     E = triu(E) + triu(E,1)';
%     
%     E=E'*E;
%     E=E-diag(diag(E))+eye(n);
%     G = E;
    %% I_e,J_e
    %%%% for fixed  diagonal entries
    I_e = [1:1:n]';
    J_e = I_e;
    k_e = length(I_e);
    
    %%  I_l,J_l
    %%%  entries with lower bounds
    I_l = [];
    J_l = [];
    for i = 1:n-ll
        r = rand(n-i,1);
        [r,ind] = sort(r);
        I_l = [I_l; i*ones(ll,1)];
        J_l = [J_l; i+ind(n-i-ll+1:n-i)];
    end
    for i = ((n-ll)+1):(n-1)
        I_l = [I_l; i*ones(n-i,1)];
        J_l = [J_l;[(i+1):n]'];
    end
    k_l = length(I_l);
    %%  I_u,J_u
    %%%%%  entries with upper bounds
    I_u = [];
    J_u = [];
    for i = 1:n-lu
        r = rand(n-i,1);
        [r,ind] = sort(r);
        I_u = [I_u; i*ones(lu,1)];
        J_u = [J_u; i+ind(n-i-lu+1:n-i)];
    end
    for i = ((n-lu)+1):(n-1)
        I_u = [I_u; i*ones(n-i,1)];
        J_u = [J_u;[(i+1):n]'];
    end
    k_u = length(I_u) ;
    %% to generate the bound e,l & u
    %%%%%%% e
    rhs    = ones(n,1);  % diagonal elements
    alpha0 = 1;          %通过调整alpha0的值来选择例子中（a）（b）两种情况，alpha0==1时对应情况为（a）
    rhs    = alpha0*rhs + (1-alpha0)*rand(n,1);
    e      = rhs;
    %%%%%%% l
    l = -0.10*rand( k_l,1)+0.6;
    %l = 0.50 * (2*rand(k_l,1)-ones(k_l,1));
    %l = 1.0 * (rand(k_l,1) - ones(k_l,1));
    
    %%%%%%% u
    u = 0.10*ones(k_u,1)+0.6;
    %u = 1.0*(rand(k_l,1) - ones(k_l,1));
    max_l = max(l);
    min_l = min(l);
    max_u = max(u);
    min_u = min(u);
    
end
end



%% realtypeG ii
%%%%%%%%%%%%%%%%%%%%%% n=1187 rank(G)=499
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
if realtypeG==2
    data1=xlsread([datapath,'\Real Examples\cormat2.csv']);
    data1=data1(:,2:end);G=data1;
    disp('Shenzhen Stocks Exchange(2016.9-2018.9), n=1187')
n=size(G,1);
%% constr_n为每行约束位置个数。
%当casetype=1时，constr_n为次对角线条数（带状对角宽度），e=1,l=-0.025r+0.5,u=+0.025r+0.5。
%当casetype=2时，constr_n为3*1向量，分别对应le,ll,lu; e=1,l=--0.3r+1,u=0.3r+0.1。
%当casetype=3时，只有下约束，ll=n，下约束为-0.1r+0.6

%% Case I  %带状上对角矩阵
%等式约束只在对角线上，且为1,（这个约束在这里无效）；
%不等式约束在次对角线上。上约束和下约束个数相等，等于次对角线条数（带宽）
%下约束为-0.025r+0.5，上约束为+0.025r+0.5
if casetype==1
fprintf('\n Real Test Example,TYPE= RE1 \n')
fprintf('\n 带状约束%d，下约束为-0.025r+0.5，上约束为+0.025r+0.5\n',constr_n(1))
bread_width=constr_n(1);
%% I_e,J_e
%%%% for fixed  diagonal entries
I_e=[1:1:n]';
J_e=I_e;
k_e=length(I_e);
%% I_l,J_l
%%%  entries with lower bounds
I_l = [];
J_l = [];
for i = 1:n-bread_width
    I_l=[I_l;i*ones(bread_width,1)];
    J_l=[J_l;(i+1:i+bread_width)'];
end
for i = ((n-bread_width)+1):(n-1)
    I_l = [I_l; i*ones(n-i,1)];
    J_l = [J_l;[(i+1):n]'];
end
k_l = length(I_l);

%% I_u,J_u
%%%  entries with lower bounds
I_u = [];
J_u = [];
for i = 1:n-bread_width
    I_u=[I_u;i*ones(bread_width,1)];
    J_u=[J_u;(i+1:i+bread_width)'];
end
for i = ((n-bread_width)+1):(n-1)
    I_u = [I_u; i*ones(n-i,1)];
    J_u = [J_u;((i+1):n)'];
end
k_u = length(I_u);
%% to generate the bound e,l & u
%%%%%%% e
rhs    = ones(n,1);  % diagonal elements
alpha0 = 1;         %通过调整alpha0的值来选择例子中（a）（b）两种情况，alpha0==1时对应情况为（a）
rhs    = alpha0*rhs + (1-alpha0)*rand(n,1);
e      = rhs;
%%%%%%% l
%l = -0.10*ones( k_l,1);
l = -0.025*rand(k_l,1)+0.5;
%l = 1.0 * (rand(k_l,1) - ones(k_l,1));

%%%%%%% u
u = 0.025*rand(k_u,1)+0.5;
%u = 1.0*(rand(k_l,1) - ones(k_l,1));
max_l = max(l);
min_l = min(l);
max_u = max(u);
min_u = min(u);
end
%% Case II  %随机位置约束
%约束位置随机
%下约束为-0.3r+0.1，上约束为0.3r+0.1
if casetype==2
    fprintf('\n Real Test Example,TYPE= RE2 \n')
    fprintf('\n 下约束为-0.3r+0.1，上约束为0.3r+0.1,约束个数为%d\n', constr_n(1))
    if length(constr_n) < 3, constr_n = repmat(constr_n(1),1,3); end
    lh =  constr_n(1);    % number of fixed off diagonal elements in each row
    ll =  constr_n(2);    % number of off diagonal elements of lower bounds in each row
    lu =  constr_n(3);    % number of off diagonal elements of upper bounds in each row
    fprintf('随机位置%d:%d:%d\n',[lh,ll,lu])
        
%     lh =  0;    % number of fixed off diagonal elements in each row
%     ll =  n;    % number of off diagonal elements of lower bounds in each row
%     lu =  0;    % number of off diagonal elements of upper bounds in each row
%     
    ll = min(ll,n-1);
    lu = min(lu,n-1);
    lh = min(lh,n-1);
    
%     E = 2.0*rand(n,n) - ones(n,n);
%     E = triu(E) + triu(E,1)';
%     
%     E=E'*E;
%     E=E-diag(diag(E))+eye(n);
%     G = E;
    %% I_e,J_e
    %%%% for fixed  diagonal entries
    I_e = [1:1:n]';
    J_e = I_e;
    k_e = length(I_e);
    
    %%  I_l,J_l
    %%%  entries with lower bounds
    I_l = [];
    J_l = [];
    for i = 1:n-ll
        r = rand(n-i,1);
        [r,ind] = sort(r);
        I_l = [I_l; i*ones(ll,1)];
        J_l = [J_l; i+ind(n-i-ll+1:n-i)];
    end
    for i = ((n-ll)+1):(n-1)
        I_l = [I_l; i*ones(n-i,1)];
        J_l = [J_l;[(i+1):n]'];
    end
    k_l = length(I_l);
    %%  I_u,J_u
    %%%%%  entries with upper bounds
    I_u = [];
    J_u = [];
    for i = 1:n-lu
        r = rand(n-i,1);
        [r,ind] = sort(r);
        I_u = [I_u; i*ones(lu,1)];
        J_u = [J_u; i+ind(n-i-lu+1:n-i)];
    end
    for i = ((n-lu)+1):(n-1)
        I_u = [I_u; i*ones(n-i,1)];
        J_u = [J_u;[(i+1):n]'];
    end
    k_u = length(I_u) ;
    %% to generate the bound e,l & u
    %%%%%%% e
    rhs    = ones(n,1);  % diagonal elements
    alpha0 = 1;          %通过调整alpha0的值来选择例子中（a）（b）两种情况，alpha0==1时对应情况为（a）
    rhs    = alpha0*rhs + (1-alpha0)*rand(n,1);
    e      = rhs;
    %%%%%%% l
    %l = -0.3*rand( k_l,1)+0.1;
    l = -0.3*rand( k_l,1)+0.1-0.075;
    
    %%%%%%% u
    %u = 0.3*rand(k_u,1)+0.1;
    u = 0.3*rand(k_u,1)+0.1+0.075;
    
    max_l = max(l);
    min_l = min(l);
    max_u = max(u);
    min_u = min(u);
    
end
%% Case III  %只有下约束且ll=n
%只有下约束为-0.1r+0.6
if casetype==3
    fprintf('\n Real Test Example,TYPE= RE3 \n')
    fprintf('\n 只有下约束为-0.1r+0.6， 所有位置%d \n',n)
    
%     lh =  constr_n(1);    % number of fixed off diagonal elements in each row
%     ll =  constr_n(2);    % number of off diagonal elements of lower bounds in each row
%     lu =  constr_n(3);    % number of off diagonal elements of upper bounds in each row
     
        
    lh =  0;    % number of fixed off diagonal elements in each row
    ll =  n;    % number of off diagonal elements of lower bounds in each row
    lu =  0;    % number of off diagonal elements of upper bounds in each row
    
    ll = min(ll,n-1);
    lu = min(lu,n-1);
    lh = min(lh,n-1);
    
%     E = 2.0*rand(n,n) - ones(n,n);
%     E = triu(E) + triu(E,1)';
%     
%     E=E'*E;
%     E=E-diag(diag(E))+eye(n);
%     G = E;
    %% I_e,J_e
    %%%% for fixed  diagonal entries
    I_e = [1:1:n]';
    J_e = I_e;
    k_e = length(I_e);
    
    %%  I_l,J_l
    %%%  entries with lower bounds
    I_l = [];
    J_l = [];
    for i = 1:n-ll
        r = rand(n-i,1);
        [r,ind] = sort(r);
        I_l = [I_l; i*ones(ll,1)];
        J_l = [J_l; i+ind(n-i-ll+1:n-i)];
    end
    for i = ((n-ll)+1):(n-1)
        I_l = [I_l; i*ones(n-i,1)];
        J_l = [J_l;[(i+1):n]'];
    end
    k_l = length(I_l);
    %%  I_u,J_u
    %%%%%  entries with upper bounds
    I_u = [];
    J_u = [];
    for i = 1:n-lu
        r = rand(n-i,1);
        [r,ind] = sort(r);
        I_u = [I_u; i*ones(lu,1)];
        J_u = [J_u; i+ind(n-i-lu+1:n-i)];
    end
    for i = ((n-lu)+1):(n-1)
        I_u = [I_u; i*ones(n-i,1)];
        J_u = [J_u;[(i+1):n]'];
    end
    k_u = length(I_u) ;
    %% to generate the bound e,l & u
    %%%%%%% e
    rhs    = ones(n,1);  % diagonal elements
    alpha0 = 1;          %通过调整alpha0的值来选择例子中（a）（b）两种情况，alpha0==1时对应情况为（a）
    rhs    = alpha0*rhs + (1-alpha0)*rand(n,1);
    e      = rhs;
    %%%%%%% l
    l = -0.10*rand( k_l,1)+0.6;
    %%%%%%% u
    u = 0.10*ones(k_u,1)+0.6;
    max_l = max(l);
    min_l = min(l);
    max_u = max(u);
    min_u = min(u);
    
end
end


%% 以下是输出
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConstrA.e = e; ConstrA.Ie = I_e; ConstrA.Je = J_e;
ConstrA.l = l; ConstrA.Il = I_l; ConstrA.Jl = J_l;
ConstrA.u = u; ConstrA.Iu = I_u; ConstrA.Ju = J_u;
end