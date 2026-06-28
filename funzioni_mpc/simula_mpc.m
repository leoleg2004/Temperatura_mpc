function [storia_x, storia_u, storia_costo] = simula_mpc(mpc_prob, x_iniziale, t_sim, Ad, Bd, dU_max, x_ref, u_ref)
    % SIMULA_MPC Risolve il problema quadratico e simula l'evoluzione del sistema
    
    % Offset affine per far sì che x_ref sia un vero punto di equilibrio
    c_affine = x_ref - Ad*x_ref - Bd*u_ref;
    nx = mpc_prob.nx;
    nu = mpc_prob.nu;
    N = mpc_prob.N;
    n_vars = mpc_prob.n_vars;
    
    storia_x = zeros(nx, t_sim+1); storia_x(:,1) = x_iniziale;
    storia_u = zeros(nu, t_sim);
    storia_costo = zeros(1, t_sim);
    u_previous = [0;0;0]; 
    options = optimoptions('quadprog', 'Display', 'off');
    
    % --- Precomputazione Matrici Invarianti (Miglioramento Performance) ---
    A_rate = zeros(nu*2*N, n_vars); 
    b_rate_base = zeros(nu*2*N, 1);
    for k = 1:N
        idx_u = (k-1)*nu+1:k*nu;
        A_rate((k-1)*2*nu+1:k*2*nu, idx_u) = [eye(nu); -eye(nu)];
        if k > 1
            idx_u_prev = (k-2)*nu+1:(k-1)*nu;
            A_rate((k-1)*2*nu+1:k*2*nu, idx_u_prev) = [-eye(nu); eye(nu)];
            b_rate_base((k-1)*2*nu+1:k*2*nu) = [dU_max; dU_max];
        end
    end
    A_ineq_tot = [mpc_prob.A_ineq_stat; A_rate];
    % ----------------------------------------------------------------------
    
    for t = 1:t_sim
        beq = repmat(c_affine, N, 1);
        beq(1:nx) = Ad * storia_x(:,t) + c_affine; 
        
        b_rate = b_rate_base;
        b_rate(1:2*nu) = [dU_max + u_previous; dU_max - u_previous];
        b_ineq_tot = [mpc_prob.b_ineq_stat; b_rate];
        
        [z_opt, fval, exitflag] = quadprog(mpc_prob.H, mpc_prob.f, A_ineq_tot, b_ineq_tot, mpc_prob.Aeq_base, beq, mpc_prob.lb, mpc_prob.ub, [], options);
        
        % Quadprog ignora la costante x_ref'*Q*x_ref. La calcoliamo matematicamente
        % come 0.5 * f' * H^-1 * f per ripristinare il VERO costo positivo!
        Cost_const = 0.5 * (mpc_prob.f' * (mpc_prob.H \ mpc_prob.f));
        storia_costo(t) = fval + Cost_const;
        
        if exitflag < 0
            error('Infeasible! Il punto allo step %d è fuori da X_N. Riduci leggermente la severità di x_iniziale.', t);
        end
        
        % L'animazione live è stata rimossa per massimizzare le prestazioni.
        
        u_applicata = z_opt(1:nu);
        storia_u(:, t) = u_applicata;
        
        storia_x(:, t+1) = Ad * storia_x(:,t) + Bd * u_applicata + c_affine;
        u_previous = u_applicata;
    end
    
    % Plot della traiettoria completa tutto in una volta alla fine
    figure('Name', 'Traiettoria 3D MPC (Temperature)', 'Color', 'w'); hold on; grid on; view(3);
    xlabel('T_1 [K]', 'Interpreter', 'latex'); 
    ylabel('T_2 [K]', 'Interpreter', 'latex'); 
    zlabel('T_3 [K]', 'Interpreter', 'latex');
    title('Evoluzione Traiettoria Temperature MPC', 'Interpreter', 'latex');
    
    plot3(storia_x(1,:), storia_x(2,:), storia_x(3,:), '-b', 'LineWidth', 2.5, 'DisplayName', 'Traiettoria Effettiva');
    
    % Plot punto di partenza
    plot3(x_iniziale(1), x_iniziale(2), x_iniziale(3), 'k*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Partenza');
    text(x_iniziale(1), x_iniziale(2), x_iniziale(3)+0.05, ' Partenza', 'FontWeight', 'bold', 'Interpreter', 'latex');
    
    % Aggiungi un pallino rosso per ogni passo effettivamente compiuto alla fine
    plot3(storia_x(1,:), storia_x(2,:), storia_x(3,:), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Passi MPC');
    
    % Plot punto di arrivo finale
    plot3(storia_x(1,end), storia_x(2,end), storia_x(3,end), 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'Arrivo Effettivo');
    text(storia_x(1,end), storia_x(2,end), storia_x(3,end)-0.05, ' Arrivo', 'FontWeight', 'bold', 'Color', 'g');
    
    % Aggiorna la legenda includendo i nuovi plot
    legend('show', 'Location', 'best', 'Interpreter', 'latex');
end
