
%% notation： this code has been embed in CONDNUM
% compute the Marix-vector(matrix) product,
% decompose the t cloumns of the latter matrix to be 1 cloumn vector 
function Ax=Axfun(flag,x,  k_e,k_l,k_u, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n)
if isequal(flag,'dim')
    Ax=k_e+k_l+k_u;
elseif isequal(flag,'real')
    Ax=1;
elseif isequal(flag,'notransp')
    for j=1:size(x,2)
        h=x(:,j);
        Ax(:,j)=Jacobian_matrix2(k_e,k_l,k_u, h, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n);
    end
elseif isequal(flag,'transp')%暂且不管转置与否
    for j=1:size(x,2)
        h=x(:,j);
        Ax(:,j)=Jacobian_matrix2(k_e,k_l,k_u, h, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n);
    end
else
    Ax=0;
end
end