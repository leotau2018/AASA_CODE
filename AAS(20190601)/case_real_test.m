%% real example
% 
% [G,ConstrA]=case_real_test(1,3,300);
% [G,ConstrA]=case_real_test(2,2,[300,300,300]);
% In our paper, we tested two cases��
%   (1): [G,ConstrA]=case_real_test(1,3,792);
%        Shanghai Stock Exchange(2016.9-2018.9), n =792, 
%        only lower bounds on all the entries: ll=n��l = -0.1r+0.6
%        denoted by E7 in the paper.
%   (2): [G,ConstrA]=case_real_test(2,2,[400,400,400);
%        Shenzhen Stocks Exchange(2016.9-2018.9), n=1187,
%        400 random positions for lower and upper bounds: e=1,l=--0.3r+1,u=0.3r+0.1��
%        denoted by E8 in the paper.
%
function [G,ConstrA]=case_real_test(realtypeG,casetype,constr_n)
%% get the path
% ���þ���·��
% datapath = 'D:\wangyunlong\��������\studysource\�о���project\3.5 AASA��Accelerated Adaptive Active Set Algorithm��';
% �������·��
datapath = '..';
%% Generate G
if realtypeG==1
%%%%%%%%%%%%%%%%%%%%%% n=792
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
data1=xlsread([datapath,'\Real Examples\cormat.csv']);
G=data1;
n=size(G,1);
%% constr_nΪÿ��Լ��λ�ø�����
%��casetype=1ʱ��constr_nΪ�ζԽ�����������״�Խǿ�ȣ���e=1,l=-0.025r+0.5,u=+0.025r+0.5��
%��casetype=2ʱ��constr_nΪ3*1�������ֱ��Ӧle,ll,lu; e=1,l=-0.025r+0.7,u=0.025r+0.7��
%��casetype=3ʱ��ֻ����Լ����ll=n����Լ��Ϊ-0.1r+0.6

%% Case I  %��״�϶ԽǾ���
%��ʽԼ��ֻ�ڶԽ����ϣ���Ϊ1,�����Լ����������Ч����
%����ʽԼ���ڴζԽ����ϡ���Լ������Լ��������ȣ����ڴζԽ�������������
%��Լ��Ϊ-0.025r+0.5����Լ��Ϊ+0.025r+0.5
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
alpha0 = 1;         %ͨ������alpha0��ֵ��ѡ�������У�a����b�����������alpha0==1ʱ��Ӧ���Ϊ��a��
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
%% Case II  %���λ��Լ��
%Լ��λ�����
%��Լ��Ϊ-0.025r+0.7����Լ��Ϊ0.025r+0.7
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
    alpha0 = 1;          %ͨ������alpha0��ֵ��ѡ�������У�a����b�����������alpha0==1ʱ��Ӧ���Ϊ��a��
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
%% Case III  %ֻ����Լ����ll=n
%��Լ��Ϊ-0.1r+0.6
if casetype==3
    fprintf('\n Real Test Example,TYPE=  RE3 \n')
    fprintf('\n ֻ����Լ��-0.1r+0.6,����λ��%d \n', n)
    
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
    alpha0 = 1;          %ͨ������alpha0��ֵ��ѡ�������У�a����b�����������alpha0==1ʱ��Ӧ���Ϊ��a��
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
%% constr_nΪÿ��Լ��λ�ø�����
%��casetype=1ʱ��constr_nΪ�ζԽ�����������״�Խǿ�ȣ���e=1,l=-0.025r+0.5,u=+0.025r+0.5��
%��casetype=2ʱ��constr_nΪ3*1�������ֱ��Ӧle,ll,lu; e=1,l=--0.3r+1,u=0.3r+0.1��
%��casetype=3ʱ��ֻ����Լ����ll=n����Լ��Ϊ-0.1r+0.6

%% Case I  %��״�϶ԽǾ���
%��ʽԼ��ֻ�ڶԽ����ϣ���Ϊ1,�����Լ����������Ч����
%����ʽԼ���ڴζԽ����ϡ���Լ������Լ��������ȣ����ڴζԽ�������������
%��Լ��Ϊ-0.025r+0.5����Լ��Ϊ+0.025r+0.5
if casetype==1
fprintf('\n Real Test Example,TYPE= RE1 \n')
fprintf('\n ��״Լ��%d����Լ��Ϊ-0.025r+0.5����Լ��Ϊ+0.025r+0.5\n',constr_n(1))
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
alpha0 = 1;         %ͨ������alpha0��ֵ��ѡ�������У�a����b�����������alpha0==1ʱ��Ӧ���Ϊ��a��
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
%% Case II  %���λ��Լ��
%Լ��λ�����
%��Լ��Ϊ-0.3r+0.1����Լ��Ϊ0.3r+0.1
if casetype==2
    fprintf('\n Real Test Example,TYPE= RE2 \n')
    fprintf('\n ��Լ��Ϊ-0.3r+0.1����Լ��Ϊ0.3r+0.1,Լ������Ϊ%d\n', constr_n(1))
    if length(constr_n) < 3, constr_n = repmat(constr_n(1),1,3); end
    lh =  constr_n(1);    % number of fixed off diagonal elements in each row
    ll =  constr_n(2);    % number of off diagonal elements of lower bounds in each row
    lu =  constr_n(3);    % number of off diagonal elements of upper bounds in each row
    fprintf('���λ��%d:%d:%d\n',[lh,ll,lu])
        
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
    alpha0 = 1;          %ͨ������alpha0��ֵ��ѡ�������У�a����b�����������alpha0==1ʱ��Ӧ���Ϊ��a��
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
%% Case III  %ֻ����Լ����ll=n
%ֻ����Լ��Ϊ-0.1r+0.6
if casetype==3
    fprintf('\n Real Test Example,TYPE= RE3 \n')
    fprintf('\n ֻ����Լ��Ϊ-0.1r+0.6�� ����λ��%d \n',n)
    
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
    alpha0 = 1;          %ͨ������alpha0��ֵ��ѡ�������У�a����b�����������alpha0==1ʱ��Ӧ���Ϊ��a��
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


%% ���������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConstrA.e = e; ConstrA.Ie = I_e; ConstrA.Je = J_e;
ConstrA.l = l; ConstrA.Il = I_l; ConstrA.Jl = J_l;
ConstrA.u = u; ConstrA.Iu = I_u; ConstrA.Ju = J_u;
end