%This program start from V'(x)=0 and V(x)=k1/lambda, find the value of L
%finally
N = 5;
a = [0.5, 3, 8, 13, 16];
b = [1, -0.1, -4, -8, -10.95];
lambda = 1;
sigma = 5;
k1 = 1.2;
finalpoint = 0.9999;
f_prev_prime_val0 = 0;
f_initial = k1/lambda;
y_values = linspace(0, 1, 10000); 

L_low = y_values(10);
L_high = y_values(end-10);
while true
        L_mid = (L_low + L_high) / 2;
        [~,Index] = min(abs(y_values-L_mid));
        diff_mid = computing(L_mid,N,a,b,lambda,sigma, f_initial, finalpoint,f_prev_prime_val0,y_values);
        if abs(diff_mid) < tolerance
            break;
        end
        
        if diff_mid < 0
            L_high = L_mid;
        else
            L_low = L_mid;
        end
end

L = L_mid;
disp(L);
pl(a,b,N,L,y_values,lambda,sigma,k1,finalpoint,f_initial)


function diff_mid = computing(L,N,a,b,lambda,sigma, f_initial, finalpoint,f_prev_prime_val0,y_values)
L_max_all = cell(N, 1);

Li_j = @(y, i, j) (a(i) * a(j) * y / (lambda * (a(i) + a(j)))) + ...
        ((a(j)^2 * b(i) - a(i)^2 * b(j)) / (lambda * (a(j)^2 - a(i)^2)));
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


[~,Index] = min(abs(y_values-L));

region = 1;
for m = N:-1:2
    L_max_m_val = L_max_all{m}(Index);
    if f_initial > L_max_m_val && L_max_m_val ~= -inf
        region = m;
        break
    end
end 

j_1 = region;


q = a(j_1)^2 / sigma^2;
k = (sqrt(1 + 8 * lambda ./ q) - 1) / 2;
D_1_1 = -((1 - L) / L)^k * (lambda * f_prev_prime_val0 * L^2 + ((-f_prev_prime_val0 - f_initial) * lambda + a(j_1) * k + a(j_1) + b(j_1)) * L + k * (-lambda * f_initial + b(j_1))) / (lambda * L * (2 * k + 1));

D_2_1 = -(-lambda * f_prev_prime_val0 * L^2 + (a(j_1) * k + (f_prev_prime_val0 + f_initial) * lambda - b(j_1)) * L + (k + 1) * (-lambda * f_initial + b(j_1))) * ((1 - L) / L)^(-k) / ((1 - L) * lambda * (2 * k + 1));
Li_j = @(y, i, j) (a(i) * a(j) * y / (lambda * (a(i) + a(j)))) + ((a(j)^2 * b(i) - a(i)^2 * b(j)) / (lambda * (a(j)^2 - a(i)^2)));

intersection_found0 = false;

for idx = 1:length(y_values)-1
    y = y_values(idx);
    for j = 2:N
        L_max_j_val = L_max_all{j}(idx);
        L_max_j_val_next = L_max_all{j}(idx+1);
        if ~isinf(L_max_j_val) && ~isinf(L_max_j_val_next)
            distance_current = f_initial - L_max_j_val;
            distance_next = f_initial - L_max_j_val_next;
            if distance_current * distance_next < 0
                intersection_idx0 = idx + 1;
                intersection_found0 = true;
                break;
            end
        end
    end
    if intersection_found0
        break;
    end
end


f_1 = @(y) D_1_1 * y .* (1 ./ y - 1).^(-k) + D_2_1 * y .* (1 ./ y - 1).^(k + 1) + (a(j_1) * y + b(j_1)) / lambda;
intersection_idx = Index;
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
    f_idx = {1:intersection_idx-1};
    f_j_list = {f_1(y_values)};
    intersection_y = y_values(intersection_idx);
    intersection_val = f_1(intersection_y);
    gamma_prev = intersection_y;
   
else
    f_idx = {1:length(y_values)};
    f_j_list = {f_1(y_values)};
    gamma_prev = 1;
end


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
        if f_prev_value > L_max_m_val  && L_max_m_val ~= -inf
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
        f_idx = [f_idx,{intersection_idx_prev+1:intersection_idx}];
        f_j_list = [f_j_list,{f_i(y_values)}];
        intersection_y = y_values(intersection_idx);
        intersection_val = f_i(intersection_y);
        
    else
        f_idx = [f_idx,{intersection_idx_prev: length(y_values)}] ;
        f_j_list = [f_j_list,{f_i(y_values)}];
        intersection_y = y_values(intersection_idx);
        intersection_val = f_i(intersection_y);
        gamma_prev = intersection_y;
        f_prev = f_i;
        fi_values = f_i(y_values);
        fi_prime = gradient(fi_values, y_values);

        f_prev_prime_val = fi_prime(intersection_idx);
    
        intersection_idx_prev = intersection_idx;
        f_prev_value = f_i(intersection_idx_prev);

        break;
    end

    gamma_prev = intersection_y;
    f_prev = f_i;
    fi_values = f_i(y_values);
    fi_prime = gradient(fi_values, y_values);

    f_prev_prime_val = fi_prime(intersection_idx);
    
    intersection_idx_prev = intersection_idx;
    f_prev_value = f_prev(y_values(intersection_idx_prev));
    
    i = i + 1;
end
if i == 2 && ~intersection_found
    final_val = f_1(finalpoint);
else
    for j = 1:i
        if ismember(finalpoint*length(y_values),f_idx{j})
            final_val = f_j_list{j}(finalpoint*length(y_values));
            break
        end
    end
end
target_val = (a(N) + b(N)) / lambda;
diff_mid= final_val - target_val;
end


function pl(a,b,N,L,y_values,lambda,sigma,k1,finalpoint,f_initial)
f_prev_prime_val = 0;
figure;
hold on;

L_max_all = cell(N, 1);
L_max_colors = {'#7E2F8E', 'green', 'blue', 'red', 'red'};

Li_j = @(y, i, j) (a(i) * a(j) * y / (lambda * (a(i) + a(j)))) + ...
        ((a(j)^2 * b(i) - a(i)^2 * b(j)) / (lambda * (a(j)^2 - a(i)^2)));
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

[~,Index] = min(abs(y_values-L));

region = 1;
for m = N:-1:2
    L_max_m_val = L_max_all{m}(Index);
    if f_initial > L_max_m_val && L_max_m_val ~= -inf
        region = m;
        break
    end
end 
j_1 = region;

q = a(j_1)^2 / sigma^2;
k = (sqrt(1 + 8 * lambda ./ q) - 1) / 2;
D_1_1 = -((1 - L) / L)^k * (lambda * f_prev_prime_val * L^2 + ((-f_prev_prime_val - f_initial) * lambda + a(j_1) * k + a(j_1) + b(j_1)) * L + k * (-lambda * f_initial + b(j_1))) / (lambda * L * (2 * k + 1));

D_2_1 = -(-lambda * f_prev_prime_val * L^2 + (a(j_1) * k + (f_prev_prime_val + f_initial) * lambda - b(j_1)) * L + (k + 1) * (-lambda * f_initial + b(j_1))) * ((1 - L) / L)^(-k) / ((1 - L) * lambda * (2 * k + 1));
Li_j = @(y, i, j) (a(i) * a(j) * y / (lambda * (a(i) + a(j)))) + ((a(j)^2 * b(i) - a(i)^2 * b(j)) / (lambda * (a(j)^2 - a(i)^2)));

intersection_found0 = false;

for idx = 1:length(y_values)-1
    y = y_values(idx);
    for j = 2:N
        L_max_j_val = L_max_all{j}(idx);
        L_max_j_val_next = L_max_all{j}(idx+1);
        if ~isinf(L_max_j_val) && ~isinf(L_max_j_val_next)
            distance_current = f_initial - L_max_j_val;
            distance_next = f_initial - L_max_j_val_next;
            if distance_current * distance_next < 0
                intersection_idx0 = idx + 1;
                intersection_found0 = true;
                break;
            end
        end
    end
    if intersection_found0
        break;
    end
end

if intersection_found0
    intersection_y0 = y_values(intersection_idx0);
    fprintf('Intersection found at y = %f\n', intersection_y0);
else
    disp('No intersection found.');
end

f_1 = @(y) D_1_1 * y .* (1 ./ y - 1).^(-k) + D_2_1 * y .* (1 ./ y - 1).^(k + 1) + (a(j_1) * y + b(j_1)) / lambda;
intersection_idx = Index;
intersection_found = false;
for idx = intersection_idx:length(y_values)
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
    f_idx = {1:intersection_idx-1};
    f_j_list = {f_1(y_values)};
    f1_values = f_1(y_values);
    f1_values(intersection_idx:end) = -inf;
    plot(y_values(Index:intersection_idx), f_1(y_values(Index:intersection_idx)), 'Color', 'black', 'LineWidth', 2, 'DisplayName', 'f_1(y)');
    intersection_y = y_values(intersection_idx);
    intersection_val = f_1(intersection_y);
    plot(intersection_y, intersection_val, 'ko', 'MarkerSize', 8, 'DisplayName', 'Intersection');
    disp(['Intersection1 found at y = ', num2str(intersection_y), ', value = ', num2str(intersection_val)]);
    
else
    f_idx = {1:length(y_values)};
    f_j_list = {f_1(y_values)};
    gamma_prev = 1;
    disp('No intersection found.');
    plot(y_values(Index:end), f_1(y_values(Index:end)), 'Color', 'black', 'LineWidth', 2, 'DisplayName', 'f_1(y)');
    final_value = f_1(finalpoint);
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
        if f_prev_value > L_max_m_val  && L_max_m_val ~= -inf
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
        f_idx = [f_idx,{intersection_idx_prev+1:intersection_idx}];
        f_j_list = [f_j_list,{f_i(y_values)}];
        plot(y_values(intersection_idx_prev:intersection_idx), f_i(y_values(intersection_idx_prev:intersection_idx)), ...
        'Color', 'black', 'LineWidth', 2, 'DisplayName', ['f_', num2str(i), '(y)']);
        intersection_y = y_values(intersection_idx);
        intersection_val = f_i(intersection_y);
        plot(intersection_y, intersection_val, 'ko', 'MarkerSize', 8, 'DisplayName', ['Intersection ', num2str(i)]);
        disp(['Intersection found at y = ', num2str(intersection_y), ', value = ', num2str(intersection_val)]);
    else
        f_idx = [f_idx,{intersection_idx_prev: length(y_values)}] ;
        f_j_list = [f_j_list,{f_i(y_values)}];
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
    f_prev_value = f_prev(y_values(intersection_idx_prev));
    
    i = i + 1;
end

if i == 2 && ~intersection_found
    final_val = f_1(finalpoint);
else
    for j = 1:i
        if ismember(finalpoint*length(y_values),f_idx{j})
            final_val = f_j_list{j}(finalpoint*length(y_values));

            break
        end
    end
end

target_val = (a(N) + b(N)) / lambda;
abs_diff = final_val - target_val;
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
title(['Plot of f_i(y) and L_max_i, k1 = ', num2str(k1)]);
legend('show');
grid on;
end

