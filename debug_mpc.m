startup_project();
try
    MPC; % Esegue MPC.m
catch e
    disp(e.message);
end
disp('--- DEBUGGING ---');
disp('Dimensione G_inf:');
disp(size(G_inf));
disp('Esiste un punto nel terminal set? (Verifico x_ref)');
disp(max(G_inf * x_ref - g_inf)); % Dovrebbe essere <= 0
disp('Verifico u_ref entro i limiti:');
disp([U_min, u_ref, U_max]);
disp('Max(dU_max):');
disp(dU_max);
disp('Verifico se il problema è ammissibile rilassando G_inf:');
% Test quadprog without G_inf
mpc_prob_no_term = mpc_prob;
mpc_prob_no_term.A_ineq_stat(end-size(G_inf,1)+1:end, :) = [];
mpc_prob_no_term.b_ineq_stat(end-size(g_inf,1)+1:end) = [];

nx = mpc_prob.nx;
nu = mpc_prob.nu;
N = mpc_prob.N;
n_vars = mpc_prob.n_vars;
t = 1;
c_affine = x_ref - Ad*x_ref - Bd*u_ref;
beq = repmat(c_affine, N, 1);
beq(1:nx) = Ad * x_iniziale + c_affine; 

A_rate = zeros(nu*2*N, n_vars); b_rate = zeros(nu*2*N, 1);
u_previous = [0;0;0];
for k = 1:N
    idx_u = (k-1)*nu+1:k*nu;
    A_rate((k-1)*2*nu+1:k*2*nu, idx_u) = [eye(nu); -eye(nu)];
    if k == 1
        b_rate((k-1)*2*nu+1:k*2*nu) = [dU_max + u_previous; dU_max - u_previous];
    else
        idx_u_prev = (k-2)*nu+1:(k-1)*nu;
        A_rate((k-1)*2*nu+1:k*2*nu, idx_u_prev) = [-eye(nu); eye(nu)];
        b_rate((k-1)*2*nu+1:k*2*nu) = [dU_max; dU_max];
    end
end
options = optimoptions('quadprog', 'Display', 'off');
[z_opt, ~, exitflag1] = quadprog(mpc_prob.H, mpc_prob.f, [mpc_prob.A_ineq_stat; A_rate], [mpc_prob.b_ineq_stat; b_rate], mpc_prob.Aeq_base, beq, mpc_prob.lb, mpc_prob.ub, [], options);
disp(['Exitflag CON terminal set: ', num2str(exitflag1)]);

[z_opt, ~, exitflag2] = quadprog(mpc_prob.H, mpc_prob.f, [mpc_prob_no_term.A_ineq_stat; A_rate], [mpc_prob_no_term.b_ineq_stat; b_rate], mpc_prob.Aeq_base, beq, mpc_prob.lb, mpc_prob.ub, [], options);
disp(['Exitflag SENZA terminal set: ', num2str(exitflag2)]);

disp('Verifico senza i limiti di rateo:');
[z_opt, ~, exitflag3] = quadprog(mpc_prob.H, mpc_prob.f, mpc_prob_no_term.A_ineq_stat, mpc_prob_no_term.b_ineq_stat, mpc_prob.Aeq_base, beq, mpc_prob.lb, mpc_prob.ub, [], options);
disp(['Exitflag SENZA terminal set e SENZA rate limit: ', num2str(exitflag3)]);

disp('Controlliamo lo stato x_1 predetto se u=0:');
disp(Ad * x_iniziale + c_affine);

