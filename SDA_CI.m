function [Q_new, P_new, step] = SDA_CI(A, Q, tol)

% X + A' * inv(X) * A = Q
% Q = z * I - A
% SDA success!

n = size(Q, 1);

E_old = A;
F_old = A';
Q_old = Q;
if issparse(Q)
    P_old = sparse(n, n);
else
    P_old = zeros(n);
end

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
    
    res = norm(Q_new + A' * (Q_new \ A) - Q, 1);
    res = res / (norm(Q_new, 1) + norm_A * norm_As * norm(Q_new, 1) + norm_Q);
    % res2 = norm(P_new + A * ((P_new - Q) \ A'), 1);
    % fprintf('res: %.3e, %.3e.\n', res, res2);
    if res < tol
        fprintf('SDA converged at %d-th step with residual = %.3e.\n', step, res);
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
    
    
    
    