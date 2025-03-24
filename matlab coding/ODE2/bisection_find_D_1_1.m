%This program start from x=0 and plot the graph
% Initial parameters
N = 5;
a = [0.5, 3, 8, 13, 16];
b = [1, -0.1, -4, -8, -10.95];
%a = [0.5, 3, 8, 13, 16];
%b = [1, -0.1, -4, -8, -10.95];
lambda = 1;
sigma = 5;
j_1 = 1;
D_2_1 = 0;
y_values = linspace(0, 1, 10000); 
finalpoint = 0.9999;
% Bisection method parameters
D_1_1_low = 0;
D_1_1_high = 300;
tolerance = 1e-1;
max_iterations = 100000;

% Bisection method loop
D_1_1_mid = find_optimal_D_1_1(a, b, lambda, sigma, N, D_1_1_low, D_1_1_high, tolerance, max_iterations,finalpoint);
disp(D_1_1_mid);
%plot the figure and show intersections
%pl(a,b,N,lambda,sigma,D_1_1_mid,D_2_1,y_values,finalpoint);


%find optimal D_1_1 such that satisfying the tolerance.
function D_1_1_mid = find_optimal_D_1_1(a, b, lambda, sigma, N, D_1_1_low, D_1_1_high, tolerance, max_iterations,finalpoint)
    iteration = 0;
    while iteration < max_iterations
        D_1_1_mid = (D_1_1_low + D_1_1_high) / 2;
        %fprintf('D_1_1_mid: %.40f\n', D_1_1_mid);
        [final_val, diff_mid, ~, ~] = compute_abs_diff(D_1_1_mid, a, b, lambda, sigma, N, finalpoint);
        
        if abs(diff_mid) < tolerance
            %fprintf('Found D_1_1 = %.40f\n with abs_diff = %e\n', D_1_1_mid, diff_mid);
            break;
        end
        
        if diff_mid < 0
            D_1_1_low = D_1_1_mid;
        else
            D_1_1_high = D_1_1_mid;
        end
        fprintf('%.40f\n', D_1_1_mid);
        fprintf('%.10f\n', diff_mid);
        iteration = iteration + 1;
    end
end


% Function to calculate final_val and abs_diff
function [final_val, diff, L_max_all, f_all] = compute_abs_diff(D_1_1, a, b, lambda, sigma, N,finalpoint)
    y_values = linspace(0, 1, 100000);
    D_2_1 = 0;
    
    Li_j = @(y, i, j) (a(i) * a(j) * y / (lambda * (a(i) + a(j)))) + ...
        ((a(j)^2 * b(i) - a(i)^2 * b(j)) / (lambda * (a(j)^2 - a(i)^2)));
    
    L_max_all = cell(N, 1);
    
    for j = 2:N
        L_max_j = -inf * ones(1, length(y_values));
        for i = 1:j-1
            current_line = arrayfun(@(y) Li_j(y, i, j), y_values);
            L_max_j = max(L_max_j, current_line);
        end
        L_max_all{j} = L_max_j;
    end

    for j = N-1:-1:2
        L_max_j = L_max_all{j};
        L_max_jp1 = L_max_all{j+1};
        exceed_indices = L_max_j > L_max_jp1;
        L_max_j(exceed_indices) = L_max_jp1(exceed_indices);
        L_max_all{j} = L_max_j;
    end
    for j = 2:N-1
        for i = j+1:N
            L_max_i = L_max_all{i};
            L_max_j = L_max_all{j};
            common_indices = find(ismember(L_max_i, L_max_j));
            L_max_j(common_indices) = -inf;
            L_max_all{j} = L_max_j;
        end
    end
    j_1 = 1;
    q = a(j_1)^2 / sigma^2;
    k = (sqrt(1 + 8 * lambda ./ q) - 1) / 2;

    f_1 = @(y) D_1_1 * y .* (1 ./ y - 1).^(-k) + D_2_1 * y .* (1 ./ y - 1).^(k + 1) + ...
        (a(j_1) * y + b(j_1)) / lambda;
    
    intersection_idx = 1;
    intersection_found = false;
    for idx = intersection_idx+1:length(y_values)-1
        y = y_values(idx);
        f1_val = f_1(y);
        for j = 2:N
            L_max_j_val = L_max_all{j}(idx);
            L_max_j_val_next = L_max_all{j}(idx+1);
            if ~isinf(L_max_j_val) && ~isinf(L_max_j_val_next)
                distance_current = f1_val - L_max_j_val;
                distance_next = f_1(y_values(idx+1)) - L_max_j_val_next;
                if distance_current * distance_next < 0
                    intersection_idx = idx+1;
                    intersection_found = true;
                    break;
                end
            end
        end
        if intersection_found
            break;
        end
    end
    
    gamma_prev = y_values(intersection_idx);
    f_prev = f_1;
    f1_values = f_1(y_values);
    f1_prime = gradient(f1_values, y_values);
    
    f_prev_prime_val = f1_prime(intersection_idx);
    intersection_idx_prev = intersection_idx;
    f_prev_value = f1_values(intersection_idx_prev);
    
    i = 2;
    f_all = {f_1};
    while true
        if gamma_prev >= 1
            break;
        end
        
        % Caculate region
        region = 1;
        for m = N:-1:2
            L_max_m_val = L_max_all{m}(intersection_idx_prev);
            if f_prev_value > L_max_m_val && L_max_m_val ~= -inf
                region = m;
                break
            end
        end
        
        next_j_i = region;

        q_i = (a(next_j_i)^2) / sigma^2;
        k_i = (sqrt(1 + 8 * lambda ./ q_i) - 1) / 2;
        D_1_i = -((1 - gamma_prev) / gamma_prev)^k_i * ...
            (lambda * f_prev_prime_val * gamma_prev^2 + ...
            ((-f_prev_prime_val - f_prev(gamma_prev)) * lambda + ...
            a(next_j_i) * k_i + a(next_j_i) + b(next_j_i)) * gamma_prev + ...
            k_i * (-lambda * f_prev(gamma_prev) + b(next_j_i))) / ...
            (lambda * gamma_prev * (2 * k_i + 1));
        D_2_i = -(-lambda * f_prev_prime_val * gamma_prev^2 + ...
            (a(next_j_i) * k_i + (f_prev_prime_val + f_prev(gamma_prev)) * lambda - b(next_j_i)) * gamma_prev + ...
            (k_i + 1) * (-lambda * f_prev(gamma_prev) + b(next_j_i))) * ...
            ((1 - gamma_prev) / gamma_prev)^(-k_i) / ...
            ((1 - gamma_prev) * lambda * (2 * k_i + 1));

        f_i = @(y) D_1_i * y .* (1 ./ y - 1).^(-k_i) + D_2_i * y .* (1 ./ y - 1).^(k_i + 1) + ...
            (a(next_j_i) * y + b(next_j_i)) / lambda;

        f_all{end+1} = f_i; % Store f_i for plotting

        intersection_found = false;
        for idx = intersection_idx+1:length(y_values)-1
            y = y_values(idx);
            fi_val = f_i(y);
            for j = 2:N
                L_max_j_val = L_max_all{j}(idx);
                L_max_j_val_next = L_max_all{j}(idx+1);
                if ~isinf(L_max_j_val) && ~isinf(L_max_j_val_next)
                    distance_current = fi_val - L_max_j_val;
                    distance_next = f_i(y_values(idx+1)) - L_max_j_val_next;
                    if distance_current * distance_next < 0
                        intersection_idx = idx+1;
                        intersection_found = true;
                        break;
                    end
                end
            end
            if intersection_found
                break;
            end
        end

        if ~intersection_found
            break;
        end

        gamma_prev = y_values(intersection_idx);
        f_prev = f_i;
        fi_values = f_i(y_values);
        fi_prime = gradient(fi_values, y_values);

        f_prev_prime_val = fi_prime(intersection_idx);
        intersection_idx_prev = intersection_idx;
        f_prev_value = fi_values(intersection_idx_prev);

        i = i + 1;
    end

    final_val = f_i(finalpoint);
    target_val = (a(N) + b(N)) / lambda;
    diff = final_val - target_val;
end

%%
%plot the result
function pl(a,b,N,lambda,sigma,D_1_1,D_2_1,y_values,finalpoint)

Li_j = @(y, i, j) (a(i) * a(j) * y / (lambda * (a(i) + a(j)))) + ((a(j)^2 * b(i) - a(i)^2 * b(j)) / (lambda * (a(j)^2 - a(i)^2)));

figure;
hold on;

L_max_all = cell(N, 1);
L_max_colors = {'#7E2F8E', 'green', 'blue', 'red', 'red'};

for j = 2:N
    L_max_j = -inf * ones(1, length(y_values));
    for i = 1:j-1
        current_line = arrayfun(@(y) Li_j(y, i, j), y_values);
        L_max_j = max(L_max_j, current_line);
    end
    L_max_all{j} = L_max_j;
end

for j = N-1:-1:2
    L_max_j = L_max_all{j};
    L_max_jp1 = L_max_all{j+1};
    exceed_indices = L_max_j > L_max_jp1;
    L_max_j(exceed_indices) = L_max_jp1(exceed_indices);
    L_max_all{j} = L_max_j;
end
for j = 2:N-1
    for i = j+1:N
        L_max_i = L_max_all{i};
        L_max_j = L_max_all{j};
        common_indices = find(ismember(L_max_i, L_max_j));
        L_max_j(common_indices) = -inf;
        L_max_all{j} = L_max_j;
    end
end
for j = 2:N
    L_max_j = L_max_all{j};
    plot(y_values, L_max_j, 'Color', L_max_colors{j-1}, 'LineWidth', 1.5);
end

plot(0, b(1) / lambda, 'ko', 'MarkerSize', 8, 'DisplayName', 'Initial Point');
plot(1, (a(N) + b(N)) / lambda, 'ko', 'MarkerSize', 8, 'DisplayName', 'Final Point');

j_1 = 1;

q = a(j_1)^2 / sigma^2;
k = (sqrt(1 + 8 * lambda ./ q) - 1) / 2;

f_1 = @(y) D_1_1 * y .* (1 ./ y - 1).^(-k) + D_2_1 * y .* (1 ./ y - 1).^(k + 1) + (a(j_1) * y + b(j_1)) / lambda;
intersection_idx = 1;
intersection_found = false;
for idx = intersection_idx+1:length(y_values)-1
    y = y_values(idx);
    f1_val = f_1(y);
    for j = 2:N
        L_max_j_val = L_max_all{j}(idx);
        L_max_j_val_next = L_max_all{j}(idx+1);
        if ~isinf(L_max_j_val) && ~isinf(L_max_j_val_next)
            distance_current = f1_val - L_max_j_val;
            distance_next = f_1(y_values(idx+1)) - L_max_j_val_next;
            if distance_current * distance_next < 0
                intersection_idx = idx+1;
                intersection_found = true;
                break;
            end
        end
    end
    if intersection_found
        break;
    end
end

if intersection_found
    plot(y_values(1:intersection_idx), f_1(y_values(1:intersection_idx)), 'Color', 'black', 'LineWidth', 2, 'DisplayName', 'f_1(y)');
    intersection_y = y_values(intersection_idx);
    intersection_val = f_1(intersection_y);
    plot(intersection_y, intersection_val, 'ko', 'MarkerSize', 8, 'DisplayName', 'Intersection');
    disp(['Intersection1 found at y = ', num2str(intersection_y), ', value = ', num2str(intersection_val)]);
else
    disp('No intersection found.');
end

gamma_prev = intersection_y;
f_prev = f_1;
f1_values = f_1(y_values);
f1_prime = gradient(f1_values, y_values);

f_prev_prime_val = f1_prime(intersection_idx);
intersection_idx_prev = intersection_idx;
f_prev_value = f1_values(intersection_idx_prev);

i = 2;
while true
    if gamma_prev >= 1
        break;
    end

    region = 1;
    for m = N:-1:2
        L_max_m_val = L_max_all{m}(intersection_idx_prev);
        if f_prev_value > L_max_m_val && L_max_m_val ~= -inf
            region = m;
            break
        end
    end
    
    next_j_i = region;

    q_i = (a(next_j_i)^2) / sigma^2;
    k_i = (sqrt(1 + 8 * lambda ./ q_i) - 1) / 2;

    D_1_i = -((1 - gamma_prev) / gamma_prev)^k_i * (lambda * f_prev_prime_val * gamma_prev^2 + ((-f_prev_prime_val - f_prev(gamma_prev)) * lambda + a(next_j_i) * k_i + a(next_j_i) + b(next_j_i)) * gamma_prev + k_i * (-lambda * f_prev(gamma_prev) + b(next_j_i))) / (lambda * gamma_prev * (2 * k_i + 1));

    D_2_i = -(-lambda * f_prev_prime_val * gamma_prev^2 + (a(next_j_i) * k_i + (f_prev_prime_val + f_prev(gamma_prev)) * lambda - b(next_j_i)) * gamma_prev + (k_i + 1) * (-lambda * f_prev(gamma_prev) + b(next_j_i))) * ((1 - gamma_prev) / gamma_prev)^(-k_i) / ((1 - gamma_prev) * lambda * (2 * k_i + 1));

    f_i = @(y) D_1_i * y .* (1 ./ y - 1).^(-k_i) + D_2_i * y .* (1 ./ y - 1).^(k_i + 1) + (a(next_j_i) * y + b(next_j_i)) / lambda;

    intersection_found = false;
    for idx = intersection_idx+1:length(y_values)-1
        y = y_values(idx);
        fi_val = f_i(y);
        for j = 2:N
            L_max_j_val = L_max_all{j}(idx);
            L_max_j_val_next = L_max_all{j}(idx+1);
            if ~isinf(L_max_j_val) && ~isinf(L_max_j_val_next)
                distance_current = fi_val - L_max_j_val;
                distance_next = f_i(y_values(idx+1)) - L_max_j_val_next;
                if distance_current * distance_next < 0
                    intersection_idx = idx+1;
                    intersection_found = true;
                    break;
                end
            end
        end
        if intersection_found
            break;
        end
    end

    if intersection_found
        plot(y_values(intersection_idx_prev:intersection_idx), f_i(y_values(intersection_idx_prev:intersection_idx)), ...
        'Color', 'black', 'LineWidth', 2, 'DisplayName', ['f_', num2str(i), '(y)']);
        intersection_y = y_values(intersection_idx);
        intersection_val = f_i(intersection_y);
        plot(intersection_y, intersection_val, 'ko', 'MarkerSize', 8, 'DisplayName', ['Intersection ', num2str(i)]);
        disp(['Intersection found at y = ', num2str(intersection_y), ', value = ', num2str(intersection_val)]);
    else
        disp(['No intersection found for f_', num2str(i), '.']);

        intersection_y = y_values(intersection_idx);
        intersection_val = f_i(intersection_y);
        gamma_prev = intersection_y;
        f_prev = f_i;
        fi_values = f_i(y_values);
        fi_prime = gradient(fi_values, y_values);

        f_prev_prime_val = fi_prime(intersection_idx);
    
        intersection_idx_prev = intersection_idx;
        f_prev_value = f_i(intersection_idx_prev);
        plot(y_values(intersection_idx_prev:end-10), f_i(y_values(intersection_idx_prev:end-10)), ...
        'Color', 'black', 'LineWidth', 2, 'DisplayName', ['f_', num2str(i), '(y)']);

        break;
    end

    gamma_prev = intersection_y;
    f_prev = f_i;
    fi_values = f_i(y_values);
    fi_prime = gradient(fi_values, y_values);

    f_prev_prime_val = fi_prime(intersection_idx);
    
    intersection_idx_prev = intersection_idx;
    f_prev_value = fi_values(intersection_idx_prev);
    
    i = i + 1;
end

final_val = f_i(finalpoint);
target_val = (a(N) + b(N)) / lambda;
abs_diff = abs(final_val - target_val);
fprintf('The absolute difference between the last value of f_i at y = 0.999 and the target is %e\n', abs_diff);
xlim([0 1])
ylim([-3 5])
text(0.4, 3, 'C_5', 'Color', 'red', 'FontSize', 10);
text(0.85, 3.3, 'C_4', 'Color', 'blue', 'FontSize', 10);
text(0.93, 2.8, 'C_3', 'Color', 'green', 'FontSize', 10);
text(0.8, 1.7, 'C_2', 'Color', '#7E2F8E', 'FontSize', 10);
text(0.5, 0, 'C_1', 'Color', '#A2142F', 'FontSize', 10);
hold off;
xlabel('y');
ylabel('f(y)');
title('Plot of f_i(y) and L_max_i');
legend('show');
grid on;
end
