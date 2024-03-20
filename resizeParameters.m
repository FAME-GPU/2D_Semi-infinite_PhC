function parameter = resizeParameters(parameter)

switch parameter.math.equation
    case 'TE'
        parameter.math.tensor22 = parameter.math.varepsilon(1 : 2, 1 : 2);
        parameter.math.tensor11 = parameter.math.mu(3, 3);
    case 'TM'
        parameter.math.tensor22 = parameter.math.mu(1 : 2, 1 : 2);
        parameter.math.tensor11 = parameter.math.varepsilon(3, 3);
end
parameter.math.inv_tensor22 = inv(parameter.math.tensor22);
parameter.math.inv_tensor11 = 1.0 / parameter.math.tensor11;

% domain 1 means air 
parameter.math.tensor22_fun = @(x)((parameter.math.tensor22 - eye(2)) * (x ~= 1) + eye(2));
parameter.math.tensor11_fun = @(x)((parameter.math.tensor11 - 1) * (x ~= 1) + 1);
parameter.math.inv_tensor22_fun = @(x)((parameter.math.inv_tensor22 - eye(2)) * (x ~= 1) + eye(2));
parameter.math.inv_tensor11_fun = @(x)((parameter.math.inv_tensor11 - 1) * (x ~= 1) + 1);

end


% if size(parameter.math.varepsilon, 1) == 1
%     parameter.math.varepsilon = parameter.math.varepsilon * eye(3);
% elseif size(parameter.math.varepsilon, 1) == 3
%     ppp = 1;
% else
%     error('Wrong dimension for permittivity!');
% end
% % domain 1 means air 
% parameter.math.varepsilon_xy_fun = @(x)((parameter.math.varepsilon(1 : 2, 1 : 2) - eye(2)) * (x == 2) + eye(2));
% parameter.math.varepsilon_zz_fun = @(x)((parameter.math.varepsilon(3, 3) - 1) * (x == 2) + 1);
% if size(parameter.math.mu, 1) == 1
%     parameter.math.mu = parameter.math.mu * eye(3);
% elseif size(parameter.math.mu, 1) == 3
%     ppp = 1;
% else
%     error('Wrong dimension for permeablity!');
% end
% parameter.math.mu_xy_fun = @(x)((parameter.math.mu(1 : 2, 1 : 2) - eye(2)) * (x == 2) + eye(2));
% parameter.math.mu_zz_fun = @(x)((parameter.math.mu(3, 3) - 1) * (x == 2) + 1);
% if size(parameter.math.varepsilon, 1) == 3
%     parameter.math.inv_varepsilon = inv(parameter.math.varepsilon);
% else
%     parameter.math.inv_varepsilon = 1 / parameter.math.varepsilon;
% end
% if size(parameter.math.mu, 1) == 3
%     parameter.math.inv_mu = inv(parameter.math.mu);
% else
%     parameter.math.inv_mu = 1 / parameter.math.mu;
% end