function [storia_x, storia_u] = simula_mpc(mpc_prob, x_iniziale, t_sim, Ad, Bd, dU_max, x_ref, u_ref)
    % SIMULA_MPC Risolve il problema quadratico e simula l'evoluzione del sistema
    
    % Offset affine per far sì che x_ref sia un vero punto di equilibrio
    c_affine = x_ref - Ad*x_ref - Bd*u_ref;
    nx = mpc_prob.nx;
    nu = mpc_prob.nu;
    N = mpc_prob.N;
    n_vars = mpc_prob.n_vars;
    
    storia_x = zeros(nx, t_sim+1); storia_x(:,1) = x_iniziale;
    storia_u = zeros(nu, t_sim);
    u_previous = [0;0;0]; 
    options = optimoptions('quadprog', 'Display', 'off');
    
    for t = 1:t_sim
        beq = repmat(c_affine, N, 1);
        beq(1:nx) = Ad * storia_x(:,t) + c_affine; 
        
        A_rate = zeros(nu*2*N, n_vars); b_rate = zeros(nu*2*N, 1);
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
        [z_opt, ~, exitflag] = quadprog(mpc_prob.H, mpc_prob.f, [mpc_prob.A_ineq_stat; A_rate], [mpc_prob.b_ineq_stat; b_rate], mpc_prob.Aeq_base, beq, mpc_prob.lb, mpc_prob.ub, [], options);
        
        if exitflag < 0
            error('Infeasible! Il punto allo step %d è fuori da X_N. Riduci leggermente la severità di x_iniziale.', t);
        end
        
        % Estrai la sequenza predetta
        X_pred = zeros(nx, N+1);
        X_pred(:, 1) = storia_x(:,t);
        for k = 1:N
            X_pred(:, k+1) = z_opt(N*nu + (k-1)*nx + 1 : N*nu + k*nx);
        end

        if t == 1
            % Inizializza i plot al primo passo
            figure('Name', 'Traiettoria 3D MPC (Temperature)', 'Color', 'w'); hold on; grid on; view(3);
            xlabel('T_1 [K]'); ylabel('T_2 [K]'); zlabel('T_3 [K]');
            h_pred = plot3(X_pred(1,:), X_pred(2,:), X_pred(3,:), '-rs', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', 'Predizione Ottima');
            h_traj = plot3(storia_x(1,1:t), storia_x(2,1:t), storia_x(3,1:t), '-b', 'LineWidth', 2.5, 'DisplayName', 'Traiettoria Effettiva');
            plot3(x_iniziale(1), x_iniziale(2), x_iniziale(3), 'k*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Partenza');
            text(x_iniziale(1), x_iniziale(2), x_iniziale(3)+0.05, ' Partenza', 'FontWeight', 'bold');
        else
            % Aggiorna dinamicamente i plot ai passi successivi
            set(h_pred, 'XData', X_pred(1,:), 'YData', X_pred(2,:), 'ZData', X_pred(3,:));
            set(h_traj, 'XData', storia_x(1,1:t), 'YData', storia_x(2,1:t), 'ZData', storia_x(3,1:t));
            drawnow limitrate;
        end
        
        u_applicata = z_opt(1:nu);
        storia_u(:, t) = u_applicata;
        
        storia_x(:, t+1) = Ad * storia_x(:,t) + Bd * u_applicata + c_affine;
        u_previous = u_applicata;
    end
    
    % Assicurati che l'ultima posizione della traiettoria sia aggiornata
    set(h_traj, 'XData', storia_x(1,:), 'YData', storia_x(2,:), 'ZData', storia_x(3,:));
    
    % Aggiungi un pallino rosso per ogni passo effettivamente compiuto alla fine
    plot3(storia_x(1,:), storia_x(2,:), storia_x(3,:), 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'DisplayName', 'Passi MPC');
    
    % Plot punto di arrivo finale
    plot3(storia_x(1,end), storia_x(2,end), storia_x(3,end), 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'Arrivo Effettivo');
    text(storia_x(1,end), storia_x(2,end), storia_x(3,end)-0.05, ' Arrivo', 'FontWeight', 'bold', 'Color', 'g');
    
    % Aggiorna la legenda includendo i nuovi plot
    legend('show', 'Location', 'best', 'Interpreter', 'latex');
end
