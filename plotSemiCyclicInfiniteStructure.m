% MATLAB代码来画一个无限重复的周期结构
clear
clc
close all

% 定义平行四边形的顶点
dx = 1; % x方向上的平移距离
dy = sqrt(3)/2; % y方向上的平移距离
parallelogram = [0 1 1.5 0.5 0; 0 0 dy dy 0];

% 定义椭圆的位置和尺寸
position = [sqrt(3)/2, dy/2-0.1, 2, 0.2]; % [x, y, width, height]

% 设置画图的数量
num_of_parallelograms = 9;

% 创建一个新的图形窗口
figure;
hold on; % 保持图形，以便在同一张图上绘制多个图形

% 绘制椭圆
position1 = position;
rectangle('Position', position1, 'EdgeColor', 'r', 'FaceColor', 'r');
position1 = position;
position1(1) = position1(1) + 3 * dx;
rectangle('Position', position1, 'EdgeColor', 'r', 'FaceColor', 'r');
position1(1) = position1(1) + 3 * dx;
rectangle('Position', position1, 'EdgeColor', 'r', 'FaceColor', 'r');

% 画出周期性结构
for i = 0 : num_of_parallelograms - 1

    % 平移平行四边形顶点
    shifted_parallelogram = bsxfun(@plus, parallelogram, [i * dx; 0]);
    
    % 绘制平行四边形
    plot(shifted_parallelogram(1, :), shifted_parallelogram(2, :), 'b-', 'LineWidth', 2);

end

shifted_parallelogram = bsxfun(@plus, parallelogram, [(i + 1) * dx; 0]);
plot(shifted_parallelogram(1, 1 : 2), shifted_parallelogram(2, 1 : 2), 'b-', 'LineWidth', 2)
plot(shifted_parallelogram(1, 3 : 4), shifted_parallelogram(2, 3 : 4), 'b-', 'LineWidth', 2)
% plot(shifted_parallelogram(1, 1 : 2) + dx, shifted_parallelogram(2, 1 : 2), 'b-', 'LineWidth', 2)
% plot(shifted_parallelogram(1, 3 : 4) + dx, shifted_parallelogram(2, 3 : 4), 'b-', 'LineWidth', 2)

for i = 0 : num_of_parallelograms - 1
    shifted_parallelogram = bsxfun(@plus, parallelogram, [i * dx; 0]);
    if mod(i, 3) == 2
        plot(shifted_parallelogram(1, 2 : 3), shifted_parallelogram(2, 2 : 3), 'm-', 'LineWidth', 2);
    end
end

% 标注方向a1和a2
a1 = [3; 0];
a2 = [1/2; sqrt(3)/2];
quiver(0, 0, a1(1), a1(2), 0, 'MaxHeadSize', 0.5/3, 'Color', '#EDB120', 'LineWidth', 2);
quiver(0, 0, a2(1), a2(2), 0, 'MaxHeadSize', 0.5, 'Color', 'm', 'LineWidth', 2);
text(a1(1) - 0.2, a1(2) + 0.2, '$\mathbf{a_1}$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', '#EDB120');
text(a2(1) - 0.5, a2(2), '$\mathbf{a_2}$', 'Interpreter', 'latex', 'FontSize', 14, 'Color', 'm');

% 用省略号表示无限延伸
text(num_of_parallelograms * dx + 0.6 * dx, 0.63 * dy, '...', 'FontSize', 20);
% text(num_of_parallelograms * dx + 1.6 * dx, 0.6 * dy, '...', 'FontSize', 20);

% 标注 0 的位置
text(-0.1, -0.1, '0', 'FontSize', 12);
  
% 标注左边界
% text(0.1 * dx, 0.6 * dy, '\partial\Omega_0', 'FontSize', 12);
% text(1.4 * dx, -0.2 * dy, '\partial\Omega_2', 'FontSize', 12);
% text(2.3 * dx, 1.2 * dy, '\partial\Omega_2', 'FontSize', 12);

% 标注区域
text(0.8 * dx, 0.8 * dy, '\Omega_0', 'FontSize', 12, 'FontWeight', 'normal');
text(1.8 * dx, 0.8 * dy, '\Omega_1', 'FontSize', 12, 'FontWeight', 'normal');
text(2.8 * dx, 0.8 * dy, '\Omega_2', 'FontSize', 12, 'FontWeight', 'normal');
text(3.8 * dx, 0.8 * dy, '\Omega_3', 'FontSize', 12, 'FontWeight', 'normal');
text(4.8 * dx, 0.8 * dy, '\Omega_1', 'FontSize', 12, 'FontWeight', 'normal');
text(5.8 * dx, 0.8 * dy, '\Omega_2', 'FontSize', 12, 'FontWeight', 'normal');
text(6.8 * dx, 0.8 * dy, '\Omega_3', 'FontSize', 12, 'FontWeight', 'normal');
text(7.8 * dx, 0.8 * dy, '\Omega_1', 'FontSize', 12, 'FontWeight', 'normal');
text(8.8 * dx, 0.8 * dy, '\Omega_2', 'FontSize', 12, 'FontWeight', 'normal');

text(0.7 * dx, 1.2 * dy, '$(K_0, M_0)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(1.7 * dx, 1.2 * dy, '$(K_1, M_1)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(2.7 * dx, 1.2 * dy, '$(K_2, M_2)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(3.7 * dx, 1.2 * dy, '$(K_3, M_3)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(4.7 * dx, 1.2 * dy, '$(K_2, M_2)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(5.7 * dx, 1.2 * dy, '$(K_3, M_3)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(6.7 * dx, 1.2 * dy, '$(K_1, M_1)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(7.7 * dx, 1.2 * dy, '$(K_2, M_2)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(8.7 * dx, 1.2 * dy, '$(K_3, M_3)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');

text(0.7 * dx, -0.2, '$(E_1, F_1)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(1.7 * dx, -0.2, '$(E_2, F_2)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(2.7 * dx, -0.2, '$(E_3, F_3)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(3.7 * dx, -0.2, '$(E_1, F_1)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(4.7 * dx, -0.2, '$(E_2, F_2)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(5.7 * dx, -0.2, '$(E_3, F_3)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(6.7 * dx, -0.2, '$(E_1, F_1)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(7.7 * dx, -0.2, '$(E_2, F_2)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
text(8.7 * dx, -0.2, '$(E_3, F_3)$', 'interpret', 'latex', 'FontSize', 12, 'FontWeight', 'normal');


% 设置坐标轴比例
axis equal;

% 扩展x轴的显示范围，以避免右侧平行四边形被截断
xlim([-0.2 * dx, num_of_parallelograms *dx + 3.2 * dx]);

% 移除坐标轴
axis off;

% 显示图形
hold off;