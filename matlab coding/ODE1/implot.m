% Parameters
a = [1, 2, 3, 4];
h = [3, 1.5, 2, 4];
U = 1;
M = 1; % V(U) = M

% Define different L and q values
L_values = [0.5, -0.8, -3];
q_values = [3.3577461243, 0.1340255737, -2.7913761139];

% Color vector
colors = ['r','g','b','m'];
labels = {'Control 1', 'Control 2', 'Control 3', 'Control 4'};
% Initialize plots
figure;
hold on;
t=0;
% Loop through each L and q pair
for i = 1:length(L_values)
    L = L_values(i);
    q = q_values(i);
    t = t+1;
    x_values = linspace(L, U, 100000);
    s = [(h(2)-h(1))/(a(2)-a(1)), (h(3)-h(2))/(a(3)-a(2)), (h(4)-h(3))/(a(4)-a(3))];
    j_1 = control(q, s);
    
    % Initial function
    f_1 = @(x) ( q - h(j_1)*L- h(j_1)/a(j_1))/a(j_1) + ...
               (-q + h(j_1)/a(j_1))/(a(j_1)*exp(-a(j_1)*L)) * exp(-a(j_1)*x) + ...
               (h(j_1)/a(j_1)) * x;
    f1_values = arrayfun(f_1, x_values);
    f1_prime = gradient(f1_values, x_values);
    
    pointfound = U;
    previndex = 1;
    current_control = 0;
    new_control = 0;
    index = length(x_values);
    findnewcontrol = false;
    
    for i = 1:length(x_values)
        % Check if control changes
        new_control = control(f1_prime(i), s);
        if new_control ~= j_1
            pointfound = x_values(i);
            index = i;
            current_control = new_control;
            findnewcontrol = true;
            break
        end
    end
    
    if ~findnewcontrol
        index = length(x_values);
    end
    
    fprev = f_1;
    fprimey = f1_prime(index);
    fy = fprev(x_values(index));
    y = x_values(index);
    
    % Plot the initial segment
    plot(x_values(previndex:index), f_1(x_values(previndex:index)), colors(j_1), 'LineWidth',4, 'HandleVisibility', 'off');
    
    j = 2;
    while index < length(x_values)
        previndex = index;
        
        newFunction_i = @(x) (-h(current_control)*y + a(current_control)*fy + fprimey - h(current_control)/a(current_control))/a(current_control) + ...
            (h(current_control)/a(current_control))*x + ...
            (-fprimey + h(current_control)/a(current_control))/(a(current_control)*exp(-a(current_control)*y)) *exp(-a(current_control)*x);
        
        f_values = arrayfun(newFunction_i, x_values);
        f_prime = gradient(f_values, x_values);
        
        findnewcontrol = false;
        for i = previndex:length(x_values)
            % Check if control changes
            new_control = control(f_prime(i), s);
            if new_control ~= current_control
                pointfound = x_values(i); 
                index = i;
                findnewcontrol = true;
                break
            end
        end
        
        if ~findnewcontrol
            index = length(x_values);
        end
        

        plot(x_values(previndex:index), newFunction_i(x_values(previndex:index)), colors(current_control), 'LineWidth', 4, 'HandleVisibility', 'off');
        current_control = new_control;
        
        fprev = newFunction_i;
        fprimey = f_prime(index);
        fy = fprev(x_values(index));
        y = x_values(index);
        
        j = j + 1;
    end
end


% Adding legend
for i = 1:length(colors)
    plot(NaN, NaN, colors(i), 'DisplayName', labels{i});
end
hold off;
%title(['Plot of V(x) with different L and q values']);
xlabel('x');
ylabel('V(x)');
legend('show', 'Location', 'best');

% Get the control
function j_1 = control(q, s)
if q >= s(3)
control = 4;
elseif q >= s(2)
control = 3;
elseif q >= s(1)
control = 2;
else
control = 1;
end
j_1 = control;
end