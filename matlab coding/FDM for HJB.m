%theta = [1.2, 0.3];
%c = [0.2, 0.4];
%delta = [0.5, 3];
theta = [1900, 1400,500];
c = [200, 400, 1100];
delta = [0.5, 4, 8];
lambda = 1;
sigma = 0.1;
q = delta.^2 ./ sigma^2;
count = length(theta);
t = 30; %The last number of V_hat displays want to select

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');


%%
% Discrete Points
zmin = -5000;
zmax = 10000;
M_min = 0; M_max = 600;  % M 的范围
M_step = 0.1;                     % M 的步长
G_start = 0;                      % G 的起始值
G_step = 0.01;                     % G 的步长
G_max = M_max;                      % G 的最大值

M_values = M_min:M_step:M_max;  % M 网格
G_values = G_start:G_step:G_max; % G 网格
num_M = length(M_values);
num_G = length(G_values);
optimal_action = NaN(num_G,num_M);
V = NaN(num_G, num_M);
%% % 边界条件 (G = 0)
for i = 1:num_M
    M = M_min + (i - 1) * M_step;
    
    % 计算所有控制变量对应的值
    values_1 = (theta + delta * M - c) / lambda;
    
    % 使用 max 函数找到最大值及其索引
    [V(1, i), j] = max(values_1);
    optimal_action(1,i) = j;
end
%%
for n = 2:num_G  % 遍历 G 网格
    G = G_values(n-1);  % 当前 G 值
    if G == 0
        for k = n:num_M-n+1
            V(2, k) = V(1, k) + G_step * 0.5 * (V(1, k+1) - 2 * V(1, k) + V(1, k-1)) / M_step^2;
            optimal_action(2,k) = optimal_action(1,k);
        end
    else
        m = num_M - n + 1;
        for k = n:m  % 遍历 M 网格（跳过边界点）
            M = M_min + (k - 1) * M_step;
            
            % 计算两个候选值对应的目标值数组
            values_2 = (-lambda * V(n-1, k) + theta + delta * M - c) ./ q;
                
            
            % 使用 max 函数找到最大值和对应的索引 w
            [~, w] = max(values_2);
            optimal_action(n,k) = w;
            % 差分计算
            term1 = (V(n-1, k+1) - 2 * V(n-1, k) + V(n-1, k-1)) / (M_step^2);
            term3 = (theta(w) + delta(w) * M - c(w) - lambda * V(n-1, k)) / (G^2 * q(w));
            V(n, k) = V(n-1, k) + G_step * (0.5 * term1 + term3);
        end
    end
end
[num_G, num_M] = size(V);

% Create meshgrid for M and G
[M_mesh, G_mesh] = meshgrid(M_values, G_values);

%% 1st graph for V 3D
colors = lines(count);  % Assign a color to each control variable
figure;
hold on;

% Iterate over each control variable value (1 to count)
for i = 1:count
    % Find positions where optimal_action == i
    mask = (optimal_action == i);
    
    % Extract corresponding M, G, V values
    M_subset = M_mesh(mask);
    G_subset = G_mesh(mask);
    V_subset = V(mask);
    
    % Plot 3D scatter plot
    scatter3(M_subset, G_subset, V_subset, 10, colors(i, :), 'filled', 'DisplayName', ['Control ', num2str(i)]);
end
theta_str = ['$\theta = [', num2str(theta(1)), ', ', num2str(theta(2)), ', ', num2str(theta(3)), ']$'];
c_str = ['$c = [', num2str(c(1)), ', ', num2str(c(2)), ', ', num2str(c(3)), ']$'];
delta_str = ['$\delta = [', num2str(delta(1)), ', ', num2str(delta(2)), ', ', num2str(delta(3)), ']$'];
lambda_str = ['$\lambda = ', num2str(lambda), '$'];
sigma_str = ['$\sigma = ', num2str(sigma), '$'];

title_text = {['$\Delta 𝕄 = ', num2str(M_step), '$, $\Delta 𝔾 = ', num2str(G_step), '$'], ...
              'Nonsmoothing', ...
              ['𝕄 \in [', num2str(M_min), ', ', num2str(M_max), ']$, $𝔾 \in [', num2str(G_start), ', ', num2str(G_max), ']$'], ...
              [theta_str, ', ', c_str, ', ',],...
              [delta_str, ', ',lambda_str, ', ', sigma_str]};
ax = gca;
xlim = ax.XLim; % 获取 x 轴范围
ylim = ax.YLim; % 获取 y 轴范围
zlim = ax.ZLim; % 获取 z 轴范围
ax = gca;
ax.XAxis.TickLabelInterpreter = 'latex';
ax.YAxis.TickLabelInterpreter = 'latex';
ax.ZAxis.TickLabelInterpreter = 'latex';

% 计算文本位置
xpos = xlim(1) - 0.1 * (xlim(2) - xlim(1))-20; % 将文本放置在 x 轴左侧
ypos = ylim(1); % 将文本放置在 y 轴底部
zpos = mean(zlim); % 将文本放置在 z 轴中间
% Set plot details
xlabel('𝕄', 'Interpreter', 'latex','FontSize', 12);
ylabel('𝔾', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'tex');
text(xpos, ypos, zpos, '$\overline{\overline{V}}$', ...
    'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Rotation', 90);

% 添加普通文本部分（(𝕄, 𝔾)）
text(xpos, ypos, zpos, ' (𝕄, 𝔾)', ...
    'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Rotation', 90);
%title(title_text, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
view(45, 30);
legend('show', 'Interpreter', 'latex', 'FontSize', 12);  % Display legend with control labels
hold off;

%% 2nd graph for Action 2D Triangle
figure;
hold on;

% Iterate over each control variable value (1 to count)
for i = 1:count
    % Find positions where optimal_action == i
    mask = (optimal_action == i);
    
    % Extract corresponding M and G values
    M_subset = M_mesh(mask);
    G_subset = G_mesh(mask);
    
    % Plot 2D scatter plot
    scatter(M_subset, G_subset, 20, colors(i, :), 'filled', 'DisplayName', ['Control ', num2str(i)]);
end

% Set plot details
xlabel('𝕄', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('𝔾', 'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
%title(title_text, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
grid on;
legend('show');  % Display legend with control labels
hold off;

%% Boundary Lines Construction
% Parameters
a = delta; % a = delta
b = theta - c; % b = theta - c

% Initialize storage for boundary_1_to_2 and boundary_2_to_3
boundary_1_to_2 = NaN(4, num_G); % Add 4th row for V_hat
boundary_2_to_3 = NaN(4, num_G); % Add 4th row for V_hat

% Loop for boundary_1_to_2 (action 1 to 2)
for g_idx = 1:num_G
    % Find the first column index where action changes from 1 to 2
    change_idx = find(optimal_action(g_idx, :) == 2, 1);
    
    % If found, extract M, G, and V values
    if ~isempty(change_idx)
        M_value = M_values(change_idx); % M value at change point
        G_value = G_values(g_idx);     % G value for this row
        V_value = V(g_idx, change_idx); % V value at change point
        
        % Calculate V_hat for action 1 to 2 (i=1, j=2)
        i = 1; j = 2; % Define indices for action 1 to 2
        V_hat = a(i)*a(j)/(lambda*(a(i)+a(j)))*M_value + ...
                (a(i)^2*b(j) - a(j)^2*b(i))/((a(i)+a(j))*(a(i)-a(j))*lambda);
        
        % Store results
        boundary_1_to_2(:, g_idx) = [M_value; G_value; V_value; V_hat];
    end
end

% Remove columns with NaN (if any G rows did not find a change point)
boundary_1_to_2 = boundary_1_to_2(:, ~any(isnan(boundary_1_to_2), 1));

% Loop for boundary_2_to_3 (action 2 to 3)
for g_idx = 1:num_G
    % Find the first column index where action changes from 2 to 3
    change_idx = find(optimal_action(g_idx, :) == 3, 1);
    
    % If found, extract M, G, and V values
    if ~isempty(change_idx)
        M_value = M_values(change_idx); % M value at change point
        G_value = G_values(g_idx);     % G value for this row
        V_value = V(g_idx, change_idx); % V value at change point
        
        % Calculate V_hat for action 2 to 3 (i=2, j=3)
        i = 2; j = 3; % Define indices for action 2 to 3
        V_hat = a(i)*a(j)/(lambda*(a(i)+a(j)))*M_value + ...
                (a(i)^2*b(j) - a(j)^2*b(i))/((a(i)+a(j))*(a(i)-a(j))*lambda);
        
        % Store results
        boundary_2_to_3(:, g_idx) = [M_value; G_value; V_value; V_hat];
    end
end

% Remove columns with NaN (if any G rows did not find a change point)
boundary_2_to_3 = boundary_2_to_3(:, ~any(isnan(boundary_2_to_3), 1));

%%
% 
% % 3rd Graph for V_hat
% % 随机选择 t 列用于展示
% num_columns_1_to_2 = size(boundary_1_to_2, 2);
% selected_columns_1_to_2 = randperm(num_columns_1_to_2, t);
% selected_boundary_1_to_2 = boundary_1_to_2(:, selected_columns_1_to_2);
% 
% num_columns_2_to_3 = size(boundary_2_to_3, 2);
% selected_columns_2_to_3 = randperm(num_columns_2_to_3, t);
% selected_boundary_2_to_3 = boundary_2_to_3(:, selected_columns_2_to_3);
% 
% % 对 boundary_1_to_2 的数据按 M (第 1 行) 进行升序排序
% [~, sort_idx_1_to_2] = sort(selected_boundary_1_to_2(1, :));
% sorted_boundary_1_to_2 = selected_boundary_1_to_2(:, sort_idx_1_to_2);
% 
% % 对 boundary_2_to_3 的数据按 M (第 1 行) 进行升序排序
% [~, sort_idx_2_to_3] = sort(selected_boundary_2_to_3(1, :));
% sorted_boundary_2_to_3 = selected_boundary_2_to_3(:, sort_idx_2_to_3);
% 
% % 转换为表格格式
% table_1 = array2table(sorted_boundary_1_to_2', ...
%     'VariableNames', {'M', 'G', 'V', 'V_hat'});
% 
% table_2 = array2table(sorted_boundary_2_to_3', ...
%     'VariableNames', {'M', 'G', 'V', 'V_hat'});
% 
% % 创建一个 figure 窗口
% figure;
% 
% % 在左侧创建 uitable 用于展示 boundary_1_to_2 的数据
% uitable('Data', table_1{:,:}, ... % 转换表格数据为数组
%     'ColumnName', table_1.Properties.VariableNames, ... % 设置列名
%     'RowName', {}, ... % 不显示行名
%     'Units', 'normalized', ... % 自适应窗口大小
%     'Position', [0.05, 0.1, 0.4, 0.8]); % 设置位置和大小 (左侧)
% 
% % 在右侧创建 uitable 用于展示 boundary_2_to_3 的数据
% uitable('Data', table_2{:,:}, ...
%     'ColumnName', table_2.Properties.VariableNames, ...
%     'RowName', {}, ...
%     'Units', 'normalized', ...
%     'Position', [0.55, 0.1, 0.4, 0.8]); % 设置位置和大小 (右侧)
% 
% % 添加标题
% 
% title_text_2 = {'Control 1 to 2 and Control 2 to 3 Tables', ...
%                 ['$\Delta M = ', num2str(M_step), '$, $\Delta G = ', num2str(G_step), '$'], ...
%               [ 'nonsmoothing'], ...
%               ['$M \in [', num2str(M_min), ', ', num2str(M_max), ']$, $G \in [', num2str(G_start), ', ', num2str(G_max), ']$']};
% 
% title(title_text_2, 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
