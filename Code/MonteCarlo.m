Initial_N = [762; 1; 0];
lambda1 = 0.003; lambda2 = 0.3;
Tf = 20;
Omega = 100;
 
time_list_list = {};
state_list_list = {};
for iter = 1:Omega
    N = Initial_N;
    
    simulated_time = 0;
    time_list = [];
    state_list = [];
    while simulated_time < Tf
        time_list = [time_list simulated_time];
        state_list = [state_list N];
        simulated_time = simulated_time - log(1 - rand(1)) / (N(1)*N(2)*lambda1 + N(2)*lambda2);
        Q = N(1) * N(2) * lambda1 / (N(1)*N(2)*lambda1 + N(2)*lambda2);
        if rand(1) < Q
            N(1) = N(1) - 1;
            N(2) = N(2) + 1;
        elseif N(2) > 0
            N(2) = N(2) - 1;
            N(3) = N(3) + 1;
        end
    end
    time_list_list = [time_list_list time_list];
    state_list_list = [state_list_list state_list];
end
 
collect_time = [];
for iter = 1:Omega
    collect_time = [collect_time time_list_list{iter}];
end
collect_time = sort(collect_time);
 
collect_state = zeros(3, length(collect_time));
for iter = 1:Omega
    state_list_list{iter} = interp1(time_list_list{iter}, transpose(state_list_list{iter}), collect_time, 'previous', 'extrap');
    collect_state = collect_state + transpose(state_list_list{iter} / Omega);
end
 
plot(collect_time, collect_state);