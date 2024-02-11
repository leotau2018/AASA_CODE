%% notation： this code has been embed in CONDNUM
% compute the Marix-vector(matrix) product,A{-1}x
% decompose the t cloumns of the latter matrix to be 1 cloumn vector 
function Ainvx=Ainverse_xfun(flag,g,  k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n)
if isequal(flag,'dim')
    Ainvx=size(g,1);
elseif isequal(flag,'real')
    Ainvx=1;
elseif isequal(flag,'notransp')
    for j=1:size(x,2)
        g=x(:,j);
        Ainvx(:,j)=pcg(@(h)Axfun(flag,h,  k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n),  g);
    end
elseif isequal(flag,'transp')%暂且不管转置与否
    for j=1:size(x,2)
        g=x(:,j);
        Ainvx(:,j)=pcg(@(h)Axfun(flag,h,  k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n),  g);
    end
else
    Ainvx=0;
end
end
