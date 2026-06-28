function mpc_prob = setup_mpc_uguaglianza(N, nx, nu, Ad, Bd, Q, R, U_min, U_max, Gx, gx, x_ref, u_ref)
    % SETUP_MPC_UGUAGLIANZA Prepara le matrici per quadprog con VINCOLO TERMINALE DI UGUAGLIANZA
    
    n_vars = N*nu + N*nx; 
    
    % Nel terminal equality constraint, si omette il costo terminale di Lyapunov
    % e si pone P = Q. La stabilità è garantita dal vincolo rigido x_N = x_ref.
    P = Q;
    
    R_blk = kron(eye(N), R);
    Q_blk = blkdiag(kron(eye(N-1), Q), P);
    H = 2 * blkdiag(R_blk, Q_blk);         
    
    % Termine lineare per l'inseguimento del riferimento
    f = zeros(n_vars, 1);                  
    for k = 1:N
        f((k-1)*nu + 1 : k*nu) = -2 * R * u_ref;
    end
    for k = 1:N-1
        f(N*nu + (k-1)*nx + 1 : N*nu + k*nx) = -2 * Q * x_ref;
    end
    f(N*nu + (N-1)*nx + 1 : N*nu + N*nx) = -2 * P * x_ref;
    
    Aeq_base = zeros(N*nx, n_vars);
    for k = 1:N
        Aeq_base((k-1)*nx + 1 : k*nx, (k-1)*nu + 1 : k*nu) = -Bd;
        Aeq_base((k-1)*nx + 1 : k*nx, N*nu + (k-1)*nx + 1 : N*nu + k*nx) = eye(nx);
        if k > 1
            Aeq_base((k-1)*nx + 1 : k*nx, N*nu + (k-2)*nx + 1 : N*nu + (k-1)*nx) = -Ad;
        end
    end
    
    lb = -inf(n_vars, 1); ub = inf(n_vars, 1); 
    for k = 1:N
        lb((k-1)*nu + 1 : k*nu) = U_min;
        ub((k-1)*nu + 1 : k*nu) = U_max;
    end
    
    % =====================================================================
    % VINCOLO TERMINALE DI UGUAGLIANZA: x_N == x_ref
    % Forziamo l'ultimo stato predetto a essere esattamente il target!
    % Questo sostituisce il Control Invariant Set (G_inf)
    % =====================================================================
    lb(N*nu + (N-1)*nx + 1 : N*nu + N*nx) = x_ref;
    ub(N*nu + (N-1)*nx + 1 : N*nu + N*nx) = x_ref;
    
    A_ineq_stat = []; b_ineq_stat = [];
    for j = 1:N-1
        A_t = zeros(size(Gx, 1), n_vars);
        A_t(:, N*nu + (j-1)*nx + 1 : N*nu + j*nx) = Gx;
        A_ineq_stat = [A_ineq_stat; A_t];
        b_ineq_stat = [b_ineq_stat; gx];
    end
    
    % Salvataggio in struttura
    mpc_prob.H = H;
    mpc_prob.f = f;
    mpc_prob.Aeq_base = Aeq_base;
    mpc_prob.lb = lb;
    mpc_prob.ub = ub;
    mpc_prob.A_ineq_stat = A_ineq_stat;
    mpc_prob.b_ineq_stat = b_ineq_stat;
    mpc_prob.n_vars = n_vars;
    mpc_prob.N = N;
    mpc_prob.nx = nx;
    mpc_prob.nu = nu;
end
