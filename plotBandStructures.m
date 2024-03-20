function [frequency, parameters] = plotBandStructures(result, parameters)

frequency = zeros(parameters.n_eigenvalues, parameters.reciprocal_lattice.n_wave_vec);
for i = 1 : parameters.reciprocal_lattice.n_wave_vec
    frequency(:, i) = result{i}.ew_finite;
end
frequency = sqrt(frequency) * parameters.lattice.a / (2 * pi);

geometries = fieldnames(parameters.geometry);
ngeom = length(geometries);
if ngeom == 3
    [~, parameters.selected_index_eigencurve2] = sort(abs(frequency(:, parameters.selected_index_wave_vec) - 0.575));
    parameters.selected_index_eigencurve = parameters.selected_index_eigencurve2(1);
end

hax = axes(figure);
for i = 1 : parameters.n_eigenvalues
    if i == parameters.selected_index_eigencurve
        plot(hax, 1 : parameters.reciprocal_lattice.n_wave_vec, frequency(i, :), 'r.-', 'MarkerSize', 10, 'LineWidth', 1), hold on
        plot(hax, parameters.selected_index_wave_vec, frequency(i, parameters.selected_index_wave_vec), 'ro', ... 
            'MarkerSize', 10, 'MarkerFaceColor', 'r'), hold on
    else
        plot(hax, 1 : parameters.reciprocal_lattice.n_wave_vec, frequency(i, :), 'b.-', 'MarkerSize', 10, 'LineWidth', 1), hold on
    end
end

if ngeom == 1
    axis([1 parameters.reciprocal_lattice.n_wave_vec 0 max(frequency(:, 1))]);
    set(hax, 'XTick', linspace(1, parameters.reciprocal_lattice.n_wave_vec, 5));
set(hax, 'XTickLabel', [-0.5, -0.25, 0, 0.25, 0.5]);
else
    axis([1 parameters.reciprocal_lattice.n_wave_vec 0.48 0.64]);
    set(hax, 'XTick', linspace(1, parameters.reciprocal_lattice.n_wave_vec, 5));
    set(hax, 'XTickLabel', [-0.2, -0.1, 0, 0.1, 0.2]);
end
xlabel('$\mathbf{k}a/(2\pi)$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$\omega a/(2\pi c)$', 'Interpreter', 'latex', 'FontSize', 24);
grid on

end