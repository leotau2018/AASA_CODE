%% Case_test����,generate the synthetic approximated matrix G and Constraints
%����:������������͡�ά����Լ������
%���:�������G��������Լ�������Ľṹ����ConstrA
function [G,ConstrA]=Case_test(casetype, n, constr_n)
%% casetype=1��2��1Ϊ��״�ԽǾ���Լ����2ΪԼ��λ������ֲ���
%% nΪ�����ά����ȡ500��1000��1500��2000
%% constr_nΪÿ��Լ��λ�ø�����
%��casetype=1ʱ��constr_nΪ�ζԽ�����������״�Խǿ�ȣ���e=1,l=-0.1,u=0.1��
%��casetype=2ʱ��constr_nΪ3*1�������ֱ��Ӧԭ����ÿһ�е�ʽ����Լ������Լ���ĸ�����e=1,l=-0.1,u=0.1��
%��casetype=3ʱ��constr_nΪ�ζԽ�����������״�Խǿ�ȣ���e=1,l=-0.2,u=0.2��
%��casetype=4ʱ��constr_nΪ3*1�������ֱ��Ӧԭ����ÿһ�е�ʽ����Լ������Լ���ĸ�����e=1,l=-0.2,u=0.2��
%��casetype=5ʱ��constr_nΪ�ζԽ�����������״�Խǿ�ȣ���e=1,l=-0.5,u=0.5��
%��casetype=6ʱ��constr_nΪ3*1�������ֱ��Ӧԭ����ÿһ�е�ʽ����Լ������Լ���ĸ�����e=1,l=-0.5,u=0.5��
%��casetype=7ʱ��constr_nΪ�ζԽ�����������״�Խǿ�ȣ���e=1,l=-0.5*rand(k_l,1),u=0.5*rang(k_u,1)��
%��casetype=8ʱ��constr_nΪ3*1�������ֱ��Ӧԭ����ÿһ�е�ʽ����Լ������Լ���ĸ�����e=1,ll=-0.5*rand(k_l,1),u=0.5*rang(k_u,1)��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if casetype==1
    disp('Test Example E1')
    constr_level = 0.1; is_random_level = 0; %  is_random_level = 0 means lower and upper bond is unrand, 1 means randomly generating
    [G, ConstrA] = diagonal_type(constr_level, n, constr_n, is_random_level);
end
if casetype==2
    disp('Test Example E2')
    constr_level = 0.1; is_random_level = 0;
    [G, ConstrA] = random_type(constr_level, n, constr_n, is_random_level);
end
if casetype==3
    disp('Test Example E3')
    constr_level = 0.2; is_random_level = 0;
    [G, ConstrA] = diagonal_type(constr_level, n, constr_n, is_random_level);
end
if casetype==4
    disp('Test Example E4')
    constr_level = 0.2; is_random_level = 0;
    [G, ConstrA] = random_type(constr_level, n, constr_n, is_random_level);
end
if casetype==5
    disp('Test Example E5')
    constr_level = 0.5; is_random_level = 0;
    [G, ConstrA] = diagonal_type(constr_level, n, constr_n, is_random_level);
end
if casetype==6
    disp('Test Example E6')
    constr_level = 0.5; is_random_level = 0;
    [G, ConstrA] = random_type(constr_level, n, constr_n, is_random_level);
end
if casetype==7
    disp('Test Example E7')
    constr_level = 0.5; is_random_level = 1;
    [G, ConstrA] = diagonal_type(constr_level, n, constr_n, is_random_level);
end
if casetype==8
    disp('Test Example E8')
    constr_level = 0.5; is_random_level = 1;
    [G, ConstrA] = random_type(constr_level, n, constr_n, is_random_level);
end
% global thetapar
% thetapar.G=G;
% thetapar.e=e;
% thetapar.I_e=I_e;
% thetapar.J_e=J_e;
% thetapar.I_l=I_l;
% thetapar.J_l=J_l;
% thetapar.I_u=I_u;
% thetapar.J_u=J_u;
% thetapar.k_e=k_e;
% thetapar.k_l=k_l;
% thetapar.k_u=k_u;
% thetapar.l=l;
% thetapar.u=u;
% thetapar.n=n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%
    function [G, ConstrA] = diagonal_type(constr_level, n, constr_n, is_random_level)
       %% ��ʽԼ��ֻ�ڶԽ����ϣ���Ϊalpha���������
        %ͨ������alpha0��ֵ��ѡ�������У�a����b�����������alpha0==1ʱ��Ӧ���Ϊ��a��
        bread_width=constr_n(1);                 %bread_widthΪ�������ԽǾ�������Խ�������(����)�ĴζԽ�������
        E = 2.0*rand(n,n) - ones(n,n);
        E = triu(E) + triu(E,1)';
        E=E'*E;
        E=E-diag(diag(E))+eye(n);
        G=E;
        % I_e,J_e
        %%% for fixed  diagonal entries
        I_e=[1:1:n]';
        J_e=I_e;
        k_e=length(I_e);
        % I_l,J_l
        %%  entries with lower bounds
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
        
        % I_u,J_u
       %%  entries with lower bounds
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
        % to generate the bound e,l & u
        %%%%%% e
        rhs    = ones(n,1);  % diagonal elements
        alpha0 = 1;         %ͨ������alpha0��ֵ��ѡ�������У�a����b�����������alpha0==1ʱ��Ӧ���Ϊ��a��
        rhs    = alpha0*rhs + (1-alpha0)*rand(n,1);
        e      = rhs;
        l = -constr_level*ones( k_l,1); %%%%%% l
        u = constr_level*ones(k_u,1); %%%%%% u
        if is_random_level  %% for random constraint level
            l = -constr_level*rand(k_l,1);
            u = constr_level*rand(k_u,1);
        end
        max_l = max(l);
        min_l = min(l);
        max_u = max(u);
        min_u = min(u);
        %% ���������
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ConstrA.e = e; ConstrA.Ie = I_e; ConstrA.Je = J_e;
        ConstrA.l = l; ConstrA.Il = I_l; ConstrA.Jl = J_l;
        ConstrA.u = u; ConstrA.Iu = I_u; ConstrA.Ju = J_u;
        % G = E;
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ���Լ��λ��
    function [G, ConstrA] = random_type(constr_level, n, constr_n, is_random_level)
        if length(constr_n)<3
            constr_n = repmat(constr_n(1),1,3);
        end
        lh =  constr_n(1);    % number of fixed off diagonal elements in each row
        ll =  constr_n(2);    % number of off diagonal elements of lower bounds in each row
        lu =  constr_n(3);    % number of off diagonal elements of upper bounds in each row
        ll = min(ll,n-1);
        lu = min(lu,n-1);
        lh = min(lh,n-1);
        
        E = 2.0*rand(n,n) - ones(n,n);
        E = triu(E) + triu(E,1)';
        
        E=E'*E;
        E=E-diag(diag(E))+eye(n);
        G = E;
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
        e      = rhs;  %%%%%% e
        l = -constr_level*ones( k_l,1); %%%%%% l
        u = constr_level*ones(k_u,1); %%%%%% u
        if is_random_level  %% for random constraint level
            l = -constr_level*rand(k_l,1);
            u = constr_level*rand(k_u,1);
        end
        max_l = max(l);
        min_l = min(l);
        max_u = max(u);
        min_u = min(u);
        %% ���������
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ConstrA.e = e; ConstrA.Ie = I_e; ConstrA.Je = J_e;
        ConstrA.l = l; ConstrA.Il = I_l; ConstrA.Jl = J_l;
        ConstrA.u = u; ConstrA.Iu = I_u; ConstrA.Ju = J_u;
        % G = E;
    end
