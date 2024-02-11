function [C,ConstrA]=gendata(n,lh,ll,lu,P)
if nargin<1
    P = 'E1';
    n=1000;    
    %%%%%%%%%%%%%%%%%%%%%%%%% Constraints
    %     lh =  0;   % number of fixed off diagonal elements in each row
    %     ll =  50;    % number of off diagonal elements of lower bounds in each row
    %     lu =  50;    % number of off diagonal elements of upper bounds in each row    
    lh=0;
    ll=200;
    lu=200;
end

switch lower(P)
    case {'e1'}
        disp('Test example E1')
        %% initiation of the problem
        %E1
        ll = min(ll,n-1);
        lu = min(lu,n-1);
        lh = min(lh,n-1);
        C=2.0*rand(n,n)-ones(n,n);
        C=triu(C)+triu(C,1)';
        for i=1:n
            C(i,i)=1;
        end
        lscalar=-0.1;uscalar=0.1;ctype=1;
        [C,ConstrA]=constrfun(C,n,lh,ll,lu,lscalar,uscalar,ctype);
        
        
    case {'e2'}
        disp('Test example E2')        
        %E2
        ll = min(ll,n-1);
        lu = min(lu,n-1);
        lh = min(lh,n-1);
        C=2.0*rand(n,n)-ones(n,n);
        C=triu(C)*triu(C,1)';
        for i=1:n
            C(i,i)=1;
        end
        lscalar=-0.1;uscalar=0.1;ctype=1;
        [C,ConstrA]=constrfun(C,n,lh,ll,lu,lscalar,uscalar,ctype);
        
        
    case {'e3'}
        disp('Test example E3')
        %E3
        ll = min(ll,n-1);
        lu = min(lu,n-1);
        lh = min(lh,n-1);
        C=2.0*rand(n,n)-ones(n,n);
        C=triu(C)+triu(C,1)';
        for i=1:n
            C(i,i)=1;
        end
        lscalar=-0.2;uscalar=0.2;ctype=2;
        [C,ConstrA]=constrfun(C,n,lh,ll,lu,lscalar,uscalar,ctype);
        
    case {'e4'}
        disp('Test example E4')
        %E4
        ll = min(ll,n-1);
        lu = min(lu,n-1);
        lh = min(lh,n-1);
        C=2.0*rand(n,n)-ones(n,n);
        C=triu(C)*triu(C,1)';
        for i=1:n
            C(i,i)=1;
        end
       lscalar=-0.2;uscalar=0.2;ctype=2;
        [C,ConstrA]=constrfun(C,n,lh,ll,lu,lscalar,uscalar,ctype);
        
    case {'e5'}
        disp('Test example E5')
        %E5
        ll = min(ll,n-1);
        lu = min(lu,n-1);
        lh = min(lh,n-1);
        C=2.0*rand(n,n)-ones(n,n);
        C=triu(C)*triu(C,1)';
        for i=1:n
            C(i,i)=1;
        end
        lscalar=0.5;uscalar=0.8;ctype=1;
        [C,ConstrA]=constrfun(C,n,lh,ll,lu,lscalar,uscalar,ctype);
        
    otherwise
        disp('Test example E6')
        %E6
        ll = min(ll,n-1);
        lu = min(lu,n-1);
        lh = min(lh,n-1);
        C=2.0*rand(n,n)-ones(n,n);
        C=triu(C)*triu(C,1)';
        for i=1:n
            C(i,i)=1;
        end
       lscalar=0.8;uscalar=0.9;ctype=2;
        [C,ConstrA]=constrfun(C,n,lh,ll,lu,lscalar,uscalar,ctype);
end



