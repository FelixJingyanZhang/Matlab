% Clinical Trial Simulation with Incremental Data for GSD

% Parameters
pi_c = 0.087; % Control group probability
risk_reduction_values = 0.15:0.01:0.22; % Percentage risk reduction (1 - rho)
rho_values = 1 - risk_reduction_values; % Relative risk (pi_e / pi_c)
n_sim = 50000; % Number of simulations
alpha = 0.025; % One-sided alpha
N_fixed = 8000; % Fixed sample size
N_max1 = 8300; % GSD1 max sample size (no futility)
N_max2 = 8750; % GSD2 max sample size (with futility)
reject_bounds = [-2.96, -2.46, -2.0]; % Rejection boundaries for GSD
futility_bounds = [0.1, -1.29, -2.0]; % Futility boundaries for GSD2
looks = [0.5, 0.7, 1.0]; % Analysis points for GSD
N_looks1 = looks * N_max1; % Sample sizes at looks for GSD1
N_looks2 = looks * N_max2; % Sample sizes at looks for GSD2

% Initialize results
power_fixed = zeros(size(rho_values));
power_gsd1 = zeros(size(rho_values));
power_gsd2 = zeros(size(rho_values));
avg_sample_fixed = zeros(size(rho_values));
avg_sample_gsd1 = zeros(size(rho_values));
avg_sample_gsd2 = zeros(size(rho_values));

% Simulation loop
for r = 1:length(rho_values)
    rho = rho_values(r);
    pi_e = pi_c * rho; % Experimental group probability
    
    % Counters for power and sample size
    reject_fixed = 0;
    reject_gsd1 = 0;
    reject_gsd2 = 0;
    sample_sizes_gsd1 = zeros(n_sim, 1);
    sample_sizes_gsd2 = zeros(n_sim, 1);
    
    for sim = 1:n_sim
        % Fixed Sample Size Design
        n_e = N_fixed / 2; % Equal allocation
        n_c = N_fixed / 2;
        y_e = binornd(n_e, pi_e); % Events in experimental group
        y_c = binornd(n_c, pi_c); % Events in control group
        pi_hat_e = y_e / n_e;
        pi_hat_c = y_c / n_c;
        if pi_hat_e == 0 || pi_hat_c == 0 || pi_hat_e == 1 || pi_hat_c == 1
            continue; % Skip to avoid log(0) or division by zero
        end
        rho_hat = pi_hat_e / pi_hat_c;
        se_ln_rho =sqrt( (2 * (1 - pi_hat_c) / (N_fixed * pi_hat_c)) + (2 * (1 - pi_hat_e) / (N_fixed*pi_hat_e)));
        z_fixed = log(rho_hat) / se_ln_rho;
        if z_fixed <= -norminv(1 - alpha)
            reject_fixed = reject_fixed + 1;
        end
        
        % GSD1 (No Futility) - Incremental Data
        stop_gsd1 = false;
        y_e_cum = 0; % Cumulative events
        y_c_cum = 0;
        for j = 1:length(looks)
            n_j = N_looks1(j);
            n_e_j = n_j / 2; % Equal allocation
            n_c_j = n_j / 2;
            % Calculate incremental sample size
            if j == 1
                delta_n_j = n_j;
            else
                delta_n_j = N_looks1(j) - N_looks1(j-1);
            end
            delta_n_e_j = delta_n_j / 2;
            delta_n_c_j = delta_n_j / 2;
            % Simulate events for incremental sample only
            y_e_new = binornd(round(delta_n_e_j), pi_e);
            y_c_new = binornd(round(delta_n_c_j), pi_c);
            % Accumulate events
            y_e_cum = y_e_cum + y_e_new;
            y_c_cum = y_c_cum + y_c_new;
            pi_hat_e_j = y_e_cum / n_e_j;
            %disp(pi_hat_e_j);
            pi_hat_c_j = y_c_cum / n_c_j;
            % if pi_hat_e_j == 0 || pi_hat_c_j == 0 || pi_hat_e_j == 1 || pi_hat_c_j == 1
            %     continue; % Skip to avoid log(0) or division by zero
            % end
            rho_hat_j = pi_hat_e_j / pi_hat_c_j;
            %disp(rho_hat_j);
            se_ln_rho_j = sqrt((2 * (1 - pi_hat_c_j) /(n_j*pi_hat_c_j) ) + (2 * (1 - pi_hat_e_j) / (n_j*pi_hat_e_j)));
            %disp(se_ln_rho_j);
            z_j = log(rho_hat_j) / se_ln_rho_j;
            if z_j <= reject_bounds(j)
                reject_gsd1 = reject_gsd1 + 1;
                sample_sizes_gsd1(sim) = n_j;
                stop_gsd1 = true;
                %disp(z_j);
                break;
                
            end
        end
        if ~stop_gsd1
            sample_sizes_gsd1(sim) = N_max1;
        end
        
        % GSD2 (With Futility) - Incremental Data
        stop_gsd2 = false;
        y_e_cum = 0; % Cumulative events
        y_c_cum = 0;
        for j = 1:length(looks)
            n_j = N_looks2(j);
            n_e_j = n_j / 2; % Equal allocation
            n_c_j = n_j / 2;
            % Calculate incremental sample size
            if j == 1
                delta_n_j = n_j;
            else
                delta_n_j = N_looks2(j) - N_looks2(j-1);
            end
            delta_n_e_j = delta_n_j / 2;
            delta_n_c_j = delta_n_j / 2;
            % Simulate events for incremental sample only
            y_e_new = binornd(round(delta_n_e_j), pi_e);
            y_c_new = binornd(round(delta_n_c_j), pi_c);
            % Accumulate events
            y_e_cum = y_e_cum + y_e_new;
            y_c_cum = y_c_cum + y_c_new;
            pi_hat_e_j = y_e_cum / n_e_j;
            pi_hat_c_j = y_c_cum / n_c_j;
            % if pi_hat_e_j == 0 || pi_hat_c_j == 0 || pi_hat_e_j == 1 || pi_hat_c_j == 1
            %     continue; % Skip to avoid log(0) or division by zero
            % end
            rho_hat_j = pi_hat_e_j / pi_hat_c_j;
            se_ln_rho_j = sqrt((2 * (1 - pi_hat_c_j) /(n_j * pi_hat_c_j) ) + (2 * (1 - pi_hat_e_j) / (n_j * pi_hat_e_j)));
            z_j = log(rho_hat_j) / se_ln_rho_j;
            if z_j <= reject_bounds(j)
                reject_gsd2 = reject_gsd2 + 1;
                sample_sizes_gsd2(sim) = n_j;
                stop_gsd2 = true;
                break;
            elseif z_j > futility_bounds(j)
                sample_sizes_gsd2(sim) = n_j;
                stop_gsd2 = true;
                break;
            end
        end
        if ~stop_gsd2
            sample_sizes_gsd2(sim) = N_max2;
        end
    end
    
    % Calculate power and average sample size
    power_fixed(r) = reject_fixed / n_sim;
    power_gsd1(r) = reject_gsd1 / n_sim;
    power_gsd2(r) = reject_gsd2 / n_sim;
    avg_sample_fixed(r) = N_fixed; % Fixed sample size is constant
    avg_sample_gsd1(r) = mean(sample_sizes_gsd1);
    avg_sample_gsd2(r) = mean(sample_sizes_gsd2);
end
%%
% Plotting Power vs. Percentage Risk Reduction
figure(1);
clf
plot(risk_reduction_values * 100, power_fixed, 'b-', 'LineWidth', 2, 'DisplayName', 'Fixed Sample');
hold on;
plot(risk_reduction_values * 100, power_gsd1, 'r-', 'LineWidth', 2, 'DisplayName', 'GSD1 (No Futility)');
plot(risk_reduction_values * 100, power_gsd2, 'g-', 'LineWidth', 2, 'DisplayName', 'GSD2 (With Futility)');
xlabel('Percentage Risk Reduction (%)');
ylabel('Power');
title({'Power vs. Percentage Risk Reduction', ...
    ['$\alpha = ', num2str(alpha), '$, $\pi_c = ', num2str(pi_c), '$, $N_{fixed} = ', num2str(N_fixed), ...
    '$, $N_{max1} = ', num2str(N_max1), '$, $N_{max2} = ', num2str(N_max2), '$']}, ...
    'Interpreter', 'latex');
legend('show', 'Location', 'best');
grid on;

% Plotting Average Sample Size vs. Percentage Risk Reduction
figure(2);
clf
plot(risk_reduction_values * 100, avg_sample_fixed, 'b-', 'LineWidth', 2, 'DisplayName', 'Fixed Sample');
hold on;
plot(risk_reduction_values * 100, avg_sample_gsd1, 'r-', 'LineWidth', 2, 'DisplayName', 'GSD1 (No Futility)');
plot(risk_reduction_values * 100, avg_sample_gsd2, 'g-', 'LineWidth', 2, 'DisplayName', 'GSD2 (With Futility)');
xlabel('Percentage Risk Reduction (%)');
ylabel('Average Sample Size');
title({'Average Sample Size vs. Percentage Risk Reduction', ...
    ['$\alpha = ', num2str(alpha), '$, $\pi_c = ', num2str(pi_c), '$, $N_{fixed} = ', num2str(N_fixed), ...
    '$, $N_{max1} = ', num2str(N_max1), '$, $N_{max2} = ', num2str(N_max2), '$']}, ...
    'Interpreter', 'latex');
legend('show', 'Location', 'best');
grid on;