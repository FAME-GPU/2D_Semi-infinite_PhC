function [E_new, F_new, Q_new, P_new, step] = SDA_cyclic_CI(A, Bs, Q, P, tol)

% X + A' * inv(X) * A = Q
% Q = z * I - A
% SDA success!

n = size(Q, 1);

E_old = A;
F_old = Bs;
Q_old = Q;
P_old = P;

step = 0;
norm_A = norm(A, 1);
norm_As = norm(A, inf);
norm_Q = norm(Q, 1);

while 1 & step < 100 
    
    step = step + 1;

    QmP = Q_old - P_old;
    E_new = E_old * (QmP \ E_old);
	F_new = F_old * (QmP \ F_old);
    Q_new = Q_old - F_old * (QmP \ E_old);
    P_new = P_old + E_old * (QmP \ F_old);

    err = norm(Q_new - Q_old, 1);
    if err < tol
        fprintf('cSDA converged at %d-th step with error = %.3e.\n', step, err);
        break;
    end

    E_old = E_new;
    F_old = F_new;
    Q_old = Q_new;
    P_old = P_new;
    
end

% if k == 100
%     fprintf('SDA did NOT converge with residual being %.3e.', res);
%     error('\n');
% end

end
    
    
    
    