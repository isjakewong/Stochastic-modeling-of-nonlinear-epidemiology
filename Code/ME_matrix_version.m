[t, y] = ode45(@mean, [0 25], [762;1;0;0;0;0;0;0;0;0;0;0]);
n1 = y(:, 1:3) + y(:, 4:6);
std1 = y(:, 1:3) + y(:, 4:6) + sqrt(abs(y(:, 7:9)));
std2 = y(:, 1:3) + y(:, 4:6) - sqrt(abs(y(:, 7:9)));

% plot temporal profile of the susceptibles population
% with their standard deviation envelop by master-equation.
plot(t,n1(:,1),'-',t,std1(:,1),'-',t,std2(:,1),'-');
title('Temporal profile of the susceptibles population.');
xlabel('Time (days)');
ylabel('Susceptibles (S)');
legend('ME mean','mean+std','mean-std')

% plot temporal profile of the infected population
% with their standard deviation envelop by master-equation.
plot(t,n1(:,2),'-',t,std1(:,2),'-',t,std2(:,2),'-');
title('Temporal profile of the infected population.');
xlabel('Time (days)');
ylabel('Infectives (I)');
legend('ME mean','mean+std','mean-std')

% plot temporal profile of the recovered population
% with their standard deviation envelop by master-equation.
plot(t,n1(:,3),'-',t,std1(:,3),'-',t,std2(:,3),'-');
title('Temporal profile of the recovered population.');
xlabel('Time (days)');
ylabel('Recovery (R)');
legend('ME mean','mean+std','mean-std')

% plot temporal profile of the S,I,R population
% without their standard deviation envelop by master-equation.
plot(t,n1(:,1),'-',t,n1(:,2),'-',t,n1(:,3),'-');
title('Temporal profile of the S,I,R population.');
xlabel('Time (days)');
ylabel('Population (N)');
legend('Susceptibles','Infectives','Recovered')

function dydt = mean(t,y)
    % function that solves the differential equations 
    % from the system size expansion of master equations
    % to simulate the temporal profile of the epidemic.
    omega = 763; lambda1 =  0.000003;
    lambda1_ = omega * lambda1; lambda2 = 0.5;
    
    phi = y(1); theta = y(2); gamma = y(3);
    e = y(4:6); v = y(7:9); c = y(10:12);
    
    dydt(1) = - lambda1_ * phi * theta;
    dydt(2) = lambda1_ * phi * theta - lambda2 * theta;
    dydt(3) = lambda2 * theta;
    
    dedt(1) = - lambda1_ * theta * e(1) - lambda1_ * phi * e(2);
    dedt(2) = lambda1_ * theta * e(1) + (lambda1_ * phi - lambda2) * e(2);
    dedt(3) = lambda2 * e(2);
   
    A = [ - theta * lambda1_, - phi * lambda1_, 0; ...,
        theta * lambda1_, (phi * lambda1_ - lambda2), 0; ...,
        0, lambda2, 0];
    B = [phi * theta * lambda1_, - phi * theta * lambda1_, 0; ...,
        - phi * theta * lambda1_, phi * theta * lambda1_ + theta * lambda2, - theta * lambda2; ...,
        0, - theta * lambda2, theta * lambda2];
    C = [v(1), c(1), c(2); c(1), v(2), c(3); c(2), c(3), v(3)];
    D = A * C + B;
    
    dydt = [dydt(1); dydt(2); dydt(3); dedt(1); dedt(2); dedt(3); ...,
        D(1,1); D(2,2); D(3,3); D(1,2) + D(2,1) - B(1,2); ...,
        D(1,3) + D(3,1) - B(1,3); D(2,3) + D(3,2) - B(2,3)];
end
