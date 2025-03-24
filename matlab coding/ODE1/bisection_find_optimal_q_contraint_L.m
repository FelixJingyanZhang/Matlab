%In this code we can set the value of L to get q, then insert it into "plotODEcontrol to get the graph"
a = [1, 2, 3, 4];
h = [3, 1.5, 2, 4];
L = -0.75;
U = 1;
M = 1; % V(U) = M
VL = 0; % V(L) = 0
tolerance = 1e-6;
colors = ['r','g','b','m'];
x_values = linspace(L, U, 1000);
s = [(h(2)-h(1))/(a(2)-a(1)), (h(3)-h(2))/(a(3)-a(2)), (h(4)-h(3))/(a(4)-a(3))]; %s = -2, 1 ,2


% Upper and lower limit manual setting.
q_lower = -10;
q_upper = 10; % V'(L) = q

% Bisection method find optimal q.
while true
    q_mid = (q_lower+q_upper)/2;
    diff = compute(a,s,h,L,U,q_mid,x_values,M);

    if abs(diff)< tolerance
        break
    end
    if diff < 0
        q_lower = q_mid;
    else
        q_upper = q_mid;
    end
end
fprintf('%.20f\n', q_mid);
disp(diff)

function diff = compute(a,s,h,L,U,q,x_values,M)
j_1 = control(q, s);
f_1 = @(x) ( q -h(j_1)*L- h(j_1)/a(j_1))/a(j_1) + ...
           (-q + h(j_1)/a(j_1))/(a(j_1)*exp(-a(j_1)*L)) * exp(-a(j_1)*x) + ...
           (h(j_1)/a(j_1)) * x;
f1_values = arrayfun(f_1, x_values);
f1_prime = gradient(f1_values, x_values);
pointfound = U;
previndex = 1;
findnewcontrol = false;

index = length(x_values);
current_control = 0;
new_control = 0;
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

    current_control = new_control;
    
    fprev = newFunction_i;
    fprimey = f_prime(index);
    fy = fprev(x_values(index));
    y = x_values(index);

    j = j + 1;

end

finalvalue = fprev(x_values(end));

diff = finalvalue - M;
end


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
