%%%%%% To generate the Jacobain product with h: theta''(z)[h] %%%%%%%
%%%%%%% 
function Ax = Jacobian_matrix2(k_e,k_l,k_u, h, I_e,J_e, I_l,J_l, I_u, J_u, Omega12,P, n)

m=k_e+k_l+k_u;
Ax =zeros(m,1);
h_e=h(1:k_e);
h_l=h(k_e+1:k_e+k_l);
h_u=h(k_e+k_l+1:m);

[r,s] = size(Omega12);
if (r>0)
    % P1 = P(:,1:r);    
    if (r< n/2)
        H = zeros(n,n);
        for i=1:k_e
            H(I_e(i), J_e(i)) = h_e(i);
        end
        for i=1:k_l
            H(I_l(i), J_l(i)) = h_l(i)+H(I_l(i), J_l(i));
        end
        for i=1:k_u
            H(I_u(i), J_u(i)) = -h_u(i)+H(I_u(i), J_u(i));  %%% upper bound
        end
        H = 0.5*(H + H');
        %H1=sparse(H)*P(:,1:r);       
        H1=H*P(:,1:r);       
        
        Omega12 = Omega12.*(H1'*P(:,r+1:n));
        
        
        %H =[(H1'*P1)*P1'+ Omega12*P2';Omega12'*P1']; %%%  H= [Omega o (P^T*diag(x)*P)]*P^T
        %H2 =[(H1'*P(:,1:r))*(P(:,1:r))'+ Omega12*(P(:,r+1:n))';Omega12'*(P(:,1:r))']; %%%  H= [Omega o (P^T*diag(x)*P)]*P^T
        %H2=P*H2;
        H = P(:,1:r)*((H1'*P(:,1:r))*(P(:,1:r))'+ 2.0*Omega12*P(:,r+1:n)');
        H=(H+H')/2;
        for i=1:k_e
            Ax(i)=H(I_e(i),J_e(i))+1e-10*h_e(i);
        end
        for i=1:k_l
            Ax(i+k_e)=H(I_l(i),J_l(i))+1e-10*h_l(i);
        end
        for i=1:k_u
            Ax(i+k_l+k_e)=-H(I_u(i),J_u(i))+1e-10*h_u(i);
        end
    else % r >n/2
        if (r==n)
            Ax =(1+1.0e-10)*h;
        else
            H = zeros(n,n);
            for i=1:k_e
                H(I_e(i), J_e(i)) = h_e(i);
            end
            for i=1:k_l
                H(I_l(i), J_l(i)) = h_l(i);
            end
            for i=1:k_u
                H(I_u(i), J_u(i)) = -h_u(i);  %%% upper bound
            end
            H = 0.5*(H + H');
            %H2=sparse(H)*P(:,r+1:n);
            H2=H*P(:,r+1:n);
            Omega12 = ones(r,s) - Omega12;
            Omega12 = Omega12.*((P(:,1:r))'*H2);
            
            %             H =[Omega12* (P(:,r+1:n))';Omega12'*(P(:,1:r))'+ ( (P(:,r+1:n))'*H2)* (P(:,r+1:n))']; %%% Assign H*P' to H= [(ones(n,n)-Omega) o (P^T*diag(x)*P)]*P^T
            %             H=P*H;
            H = P(:,r+1:n)*((H2'*P(:,r+1:n))*(P(:,r+1:n))'+ 2.0*Omega12'*P(:,1:r)');
            H=(H+H')/2;
            for i=1:k_e
                Ax(i)=h_e(i)-H(I_e(i),J_e(i))+1e-10*h_e(i);
            end
            for i=1:k_l
                Ax(i+k_e)=h_l(i)-H(I_l(i),J_l(i))+1e-10*h_l(i);
            end
            for i=1:k_u
                Ax(i+k_l+k_e)=-h_u(i)+H(I_u(i),J_u(i))+1e-10*h_u(i);
            end
        end
    end
end
return
