function [result_CI, result_select] = computeSelectedEigenpairs3(mtx, result, parameters)
% Use contour integral to compute the eigenpairs of infinite structures

K_C = mtx.K_C;
M_C = mtx.M_C;
K_L = mtx.K_L;
M_L = mtx.M_L;
E_L = mtx.E_L;
F_L = mtx.F_L;
E_LB = mtx.E_LB;
F_LB = mtx.F_LB;
K_R = mtx.K_R;
M_R = mtx.M_R;
E_R = mtx.E_R;
F_R = mtx.F_R;
E_RB = mtx.E_RB;
F_RB = mtx.F_RB;

m1 = size(K_L, 1);
m2 = size(K_C, 1);
m3 = size(K_R, 1);
num_points = parameters.contour_integral.num_points;
radius = parameters.contour_integral.radius;
tol = parameters.SDA.tol;

ew_finite_select = result{parameters.selected_index_wave_vec}.ew_finite(parameters.selected_index_eigencurve);
ev_finite_select = result{parameters.selected_index_wave_vec}.ev_finite(:, parameters.selected_index_eigencurve);

% ew_finite_select = (0.58 * 2 * pi)^2;

%% contour integral
time_contour_integral = tic;
ev = 0;
ew_ev = 0;
time_SDA = zeros(num_points, 1);
for i = 1 : num_points
    z = ew_finite_select + radius * exp(1i * 2 * pi * i / num_points);
    ACz = z * M_C - K_C;
    ALz = z * M_L - K_L;
    ARz = z * M_R - K_R;
    NLz = E_L - z * F_L;
    NRz = E_R - z * F_R;
    NLBz = E_LB - z * F_LB;
    NRBz = E_RB - z * F_RB;

    tmp_time_SDA_L = tic;
    [invGLz, ~, stepL] = SDA_CI(NLz', ALz, tol);
    [invGRz, ~, stepR] = SDA_CI(NRz, ARz, tol);
    invGCz = ACz - NLBz * (invGLz \ NLBz') - NRBz' * (invGRz \ NRBz);
    time_SDA(i) = toc(tmp_time_SDA_L);

    dt = 2 * pi / num_points;
    dz = 1i * radius * exp(1i * 2 * pi * i / num_points) * dt;
    ele_integral = invGCz \ (dz * speye(m2));
    ev = ev + ele_integral;
    ew_ev = ew_ev + z * ele_integral;
end
ev = ev / (2 * pi * 1i);
ew_ev = ew_ev / (2 * pi * 1i);
time_contour_integral = toc(time_contour_integral);

ew_one = eig(ev);
[~, index] = sort(abs(ew_one), 'descend'); 
ew_one = ew_one(index);
isRank1 = 1 - abs(ew_one(2) / ew_one(1));
fprintf('Is rank 1: %.12f.\n', isRank1);

ew_mu = eig(ew_ev);
[~, index] = sort(abs(ew_mu), 'descend'); 
ew_mu = ew_mu(index);
isRank1_2 = 1 - abs(ew_mu(2) / ew_mu(1));
fprintf('Is rank 1_2: %.12f.\n', isRank1_2);

mu = real(ew_mu(1) / ew_one(1));

[qC, ~] = qr(ev);
qC = qC(:, 1);

qC_finite = ev_finite_select(parameters.num_blocks * m1 + 1 : parameters.num_blocks * m1 + m2);
qC_finite = qC_finite / norm(qC_finite);
ratio_qC = qC_finite ./ qC;
isEigenvectorC = 1 - std(ratio_qC);
fprintf('Is eigenvector: %.12f.\n', isEigenvectorC);

result_CI.ev_inf = ev;
result_CI.ew_ev_inf = ew_ev;
result_CI.qC_inf = qC;
result_CI.ratio_qC = ratio_qC;
result_CI.isRank1 = isRank1;
result_CI.isEigenvectorC = isEigenvectorC;
result_CI.time_contour_integral = time_contour_integral;
result_CI.time_SDA = time_SDA;
result_CI.ew_one = ew_one;
result_CI.ew_mu = ew_mu;
result_CI.mu = mu;
result_CI.mu_normalized = sqrt(mu) * parameters.lattice.a / (2 * pi);

result_select.qC_finite = qC_finite;
result_select.ev_finite_select = ev_finite_select;
result_select.ew_finite_select = ew_finite_select;

end