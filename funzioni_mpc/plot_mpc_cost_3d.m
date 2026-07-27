function plot_mpc_cost_3d(mpc_prob, A_long_ds, dU_max, storia_x, storia_costo, Ts)
    % =========================================================================
    % PLOT_MPC_COST_3D - Visualizza il Funzionale di Costo Ottimo J*(x) dell'MPC
    % Adattato per il sistema delle 3 Camere
    % =========================================================================

    if nargin < 6
        Ts = 120;
    end

    t_lim = 10; % Delta K
    x_ref = storia_x(:,end);

    % Usiamo una griglia per T1 e T2 (le temperature principali)
    % IN AUTOMATION F16 GLI STATI ERANO ASSOLUTI, QUINDI CENTRIAMO LA GRIGLIA
    [T1_grid, T2_grid] = meshgrid(linspace(x_ref(1)-t_lim, x_ref(1)+t_lim, 20), linspace(x_ref(2)-t_lim, x_ref(2)+t_lim, 20));

    nx = mpc_prob.nx;
    nu = mpc_prob.nu;

    Q_lqr = blkdiag(100, 100, 100, 1, 1, 1);
    R_lqr = eye(nu) * 0.1;
    [~, P_ds, ~] = dlqr(A_long_ds, mpc_prob.Aeq_base(1:nx, 1:nu)*(-1), Q_lqr, R_lqr); 
    
    V_surf = zeros(size(T1_grid));
    for i = 1:size(T1_grid, 1)
        for j = 1:size(T1_grid, 2)
            % La conca si calcola sull'errore
            stato_surf = [T1_grid(i,j) - x_ref(1); T2_grid(i,j) - x_ref(2); 0; 0; 0; 0];
            V_surf(i,j) = stato_surf' * P_ds * stato_surf;
        end
    end
    
    max_v_traj = max(storia_costo);
    if isempty(max_v_traj) || isnan(max_v_traj)
        max_v_traj = 1;
    end
    
    % ================= SCALATURA DEL COSTO =================
    % Ripristiniamo la scala lineare come prima, dividendo per 10^4 per avere 
    % numeri leggibili sull'asse (es. 1, 2, 3 invece di 10000, 20000).
    fattore_scala = 1e4;
    V_surf = V_surf / fattore_scala;
    storia_costo = storia_costo / fattore_scala;
    max_v_traj = max_v_traj / fattore_scala;
    
    % --- RISOLUZIONE SCHIACCIAMENTO VISIVO ---
    % Calcoliamo l'altezza teorica della superficie esattamente al punto di partenza
    stato_start = [storia_x(1,1) - x_ref(1); storia_x(2,1) - x_ref(2); 0; 0; 0; 0];
    V_start_teorico = (stato_start' * P_ds * stato_start) / fattore_scala;
    
    % Scaliamo la superficie in modo che "tocchi" perfettamente la linea rossa alla partenza!
    if V_start_teorico > 0
        scale_factor = max_v_traj / V_start_teorico;
        V_surf = V_surf * scale_factor;
    end

    %% --- FIGURA 3D ---
    figure('Name', sprintf('Costo MPC 3D (Ts = %g s)', Ts), 'Color', 'w', 'Position', [150 150 850 650]);
    hold on; grid on;
    
    h_surf = surf(T1_grid, T2_grid, V_surf, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
    colormap jet;
    cb = colorbar;
    ylabel(cb, 'Costo Ottimo MPC $J^*(x_k) \times 10^4$', 'Interpreter', 'latex', 'FontSize', 12);
    
    % Traiettoria
    N_steps = length(storia_costo);
    t_discrete = 1:N_steps;
    t_fine = linspace(1, N_steps, N_steps * 10);
    
    T1_smooth = pchip(t_discrete, storia_x(1,1:N_steps), t_fine);
    T2_smooth = pchip(t_discrete, storia_x(2,1:N_steps), t_fine);
    V_smooth = pchip(t_discrete, storia_costo, t_fine);
    
    % Riduciamo anche l'offset visivo per non staccare troppo la linea dal fondo
    z_offset = max(max(V_surf(:)), max_v_traj) * 0.01;
    
    h_line = plot3(T1_smooth, T2_smooth, V_smooth + z_offset*0.5, '-r', 'LineWidth', 2);
    
    % Punti discreti veri
    plot3(storia_x(1,1:N_steps), storia_x(2,1:N_steps), storia_costo + z_offset, 'o', ...
          'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w', 'MarkerSize', 4);
          
    % Start (Punto nero preciso all'inizio del funzionale)
    h_start = plot3(storia_x(1,1), storia_x(2,1), storia_costo(1) + z_offset*2, 'ko', ...
                    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    text(storia_x(1,1), storia_x(2,1), storia_costo(1) + z_offset*4, ' Partenza', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
    
    % End
    h_end = plot3(storia_x(1,N_steps), storia_x(2,N_steps), storia_costo(N_steps) + z_offset*2, 's', ...
          'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'MarkerSize', 8);

    title(sprintf('\\textbf{Funzionale di Costo MPC $J^*(x_k)$ a Tempo Discreto ($T_s = %g$ s)}', Ts), 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('Temperatura $T_1$ [K]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Temperatura $T_2$ [K]', 'Interpreter', 'latex', 'FontSize', 12);
    zlabel('Costo $J^*(x_k) \times 10^4$', 'Interpreter', 'latex', 'FontSize', 12);
    
    legend([h_surf, h_line, h_start, h_end], ...
           {'Superficie $J^*(x_k)$', 'Evoluzione di Stato', 'Partenza ($k=0$)', 'Arrivo'}, ...
           'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
    
    zlim([0, max(max(V_surf(:)), max_v_traj) * 1.2]); 
    xlim([x_ref(1)-t_lim, x_ref(1)+t_lim]);
    ylim([x_ref(2)-t_lim, x_ref(2)+t_lim]);
    view(-35, 30);
end
