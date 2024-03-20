function [X, Y] = FPI(B, Q, c, tol)
%% Q is Hermitian
%% fix point

X = Q;
Y = inv(Q);
step = 0;

while step < 100000
    
    step = step + 1;
    
    X = c * X + (1 - c) * (Q - B' * (X \ B));
    Y = c * Y + (1 - c) * inv(Q - B' * Y * B);
    
    res = norm(X + B' * (X \ B) - Q, 1);
    res2 = norm(inv(Y) + B' * (Y * B) - Q, 1);
    % fprintf('res: %.3e, %.3e.\n', res, res2);
    if max(res, res2) < tol
        % fprintf('res: %.3e, %.3e.\n', res, res2);
        fprintf('FPI converged at the %d-th step with rate = %.2f and residual = %.3e.\n', step, c, max(res, res2));
        break;
    end
   
    if step == 10000
        fprintf('FPI did NOT converge at the %d-th step with rate = %.2f.\n', step, c);
    end
end

end



