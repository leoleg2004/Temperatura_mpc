function plot_lyapunov_discrete(P_ds, A_cl_ds, Ts, x_iniziale, x_ref)
    % =========================================================================
    % PLOT_LYAPUNOV_DISCRETE - Ritratto di Fase e Funzione di Lyapunov 3D 
    % Adattato per il progetto MPC Appartamento 3 Stanze
    % =========================================================================
    
    if nargin < 4
        x_iniziale = [284; 285; 284; 0; 10; 0];
        x_ref = [289; 289; 289; 100; 100; 100];
    end
    
    % Lavoriamo con gli errori (delta) rispetto al riferimento
    x0_delta = x_iniziale - x_ref;
    
    % Aggiungiamo anche un paio di altre condizioni iniziali per far 
    % vedere bene le curve (come in F16), più la tua vera partenza.
    x0_mult = [x0_delta, [-5; 3; 0; 0; 0; 0], [4; -6; 0; 0; 0; 0]];
    
    set(0,'DefaultLineLineWidth',1.5);
    set(0,'DefaultAxesFontSize',14);
    
    % --- 1. LIMITI VISIVI ---
    t_lim = max(abs(x0_delta(1:2))) + 2; 
    if t_lim < 10; t_lim = 10; end
    
    [T1_grid, T2_grid] = meshgrid(linspace(-t_lim, t_lim, 50), linspace(-t_lim, t_lim, 50));
    
    % --- 2. RITRATTO DI FASE 2D (Campo Vettoriale e Contour) ---
    figure('Name', 'Ritratto di Fase LQR (Discreto)', 'Color', 'w', 'Position', [50 100 800 600]);
    hold on; grid on;
    
    dT1 = zeros(size(T1_grid));
    dT2 = zeros(size(T2_grid));
    V_contour = zeros(size(T1_grid));
    
    for i = 1:numel(T1_grid)
        stato = [T1_grid(i); T2_grid(i); 0; 0; 0; 0];
        stato_next = A_cl_ds * stato;
        
        % Il vettore freccia è la differenza di stato (Delta T)
        dT1(i) = stato_next(1) - stato(1); 
        dT2(i) = stato_next(2) - stato(2); 
        V_contour(i) = stato' * P_ds * stato;
    end
    
    % Campo Vettoriale (Normalizzato per non avere frecce giganti)
    L = sqrt(dT1.^2 + dT2.^2) + 1e-6;
    h_quiv = quiver(T1_grid, T2_grid, dT1./L, dT2./L, 0.5, 'Color', [0.6 0.6 0.6]);
    
    % Curve di Livello (Spaziatura Logaritmica)
    val_min = min(V_contour(:));
    val_max = max(V_contour(:));
    livelli_log = logspace(log10(val_min + 1e-1), log10(val_max), 20);
    [~, h_cont] = contour(T1_grid, T2_grid, V_contour, livelli_log, 'LineWidth', 1.2, 'LineColor', [0.2 0.5 0.8]);
    
    % Simulazione e Plot delle Traiettorie 2D
    N_steps = 150;
    colori = lines(size(x0_mult, 2));
    
    for i = 1:size(x0_mult, 2)
        x_traj = zeros(6, N_steps);
        x_traj(:, 1) = x0_mult(:, i);
        for k = 1:N_steps-1
            x_traj(:, k+1) = A_cl_ds * x_traj(:, k);
        end
        if i == 1
            h_traj = plot(x_traj(1,:), x_traj(2,:), '-', 'Color', colori(i,:), 'LineWidth', 2.5);
            h_start = plot(x_traj(1,1), x_traj(2,1), '*', 'Color', colori(i,:), 'MarkerSize', 10, 'LineWidth', 2);
        else
            plot(x_traj(1,:), x_traj(2,:), '-', 'Color', colori(i,:), 'LineWidth', 2);
            plot(x_traj(1,1), x_traj(2,1), '*', 'Color', colori(i,:), 'MarkerSize', 8);
        end
    end
    
    xlabel('Deviazione $\Delta T_1$ [K]', 'Interpreter', 'latex') 
    ylabel('Deviazione $\Delta T_2$ [K]', 'Interpreter', 'latex') 
    title('\textbf{Ritratto di Fase LQR e Curve di Livello $V(x)$}', 'Interpreter', 'latex')
    
    xlim([-t_lim, t_lim]); 
    ylim([-t_lim, t_lim]);
    
    legend([h_traj, h_start, h_quiv, h_cont], ...
           {'Traiettoria (Vera Partenza)', 'Condizione Iniziale', 'Campo Vettoriale', 'Curve $V(x)=cost$'}, ...
           'Interpreter','latex', 'Location', 'bestoutside');


    % --- 3. FIGURA 3D DELLA FUNZIONE DI LYAPUNOV ---
    figure('Name', sprintf('Lyapunov 3D Discreto (Ts = %g s)', Ts), 'Color', 'w', 'Position', [150 150 850 650]);
    hold on; grid on;
    
    h_surf = surf(T1_grid, T2_grid, V_contour, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
    colormap jet;
    cb = colorbar;
    ylabel(cb, 'Energia $V(x_k)$', 'Interpreter', 'latex', 'FontSize', 12);
    
    max_v_traj = 0; 
    
    for i = 1:size(x0_mult, 2)
        x_traj = zeros(6, N_steps);
        x_traj(:, 1) = x0_mult(:, i);
        for k = 1:N_steps-1
            x_traj(:, k+1) = A_cl_ds * x_traj(:, k);
        end
        
        V_traj = zeros(N_steps, 1);
        for k = 1:N_steps
            V_traj(k) = x_traj(:, k)' * P_ds * x_traj(:, k);
        end
        
        if max(V_traj) > max_v_traj
            max_v_traj = max(V_traj);
        end
        
        t_discrete = 1:N_steps;
        t_fine = linspace(1, N_steps, N_steps * 10);
        T1_smooth = pchip(t_discrete, x_traj(1,:), t_fine);
        T2_smooth = pchip(t_discrete, x_traj(2,:), t_fine);
        V_smooth = pchip(t_discrete, V_traj, t_fine);
        
        if i == 1
            h_line = plot3(T1_smooth, T2_smooth, V_smooth, '-', 'Color', colori(i,:), 'LineWidth', 2.5);
        else
            plot3(T1_smooth, T2_smooth, V_smooth, '-', 'Color', colori(i,:), 'LineWidth', 2);
        end
        
        % Punti discreti
        plot3(x_traj(1,:), x_traj(2,:), V_traj, 'o', 'MarkerEdgeColor', colori(i,:), 'MarkerFaceColor', 'w', 'MarkerSize', 4);
              
        if i == 1
            h_start3 = plot3(x_traj(1,1), x_traj(2,1), V_traj(1), 'o', 'MarkerFaceColor', colori(i,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        else
            plot3(x_traj(1,1), x_traj(2,1), V_traj(1), 'o', 'MarkerFaceColor', colori(i,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        end
        
        plot3(x_traj(1,end), x_traj(2,end), V_traj(end), 's', 'MarkerFaceColor', colori(i,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    
    h_end = plot3(0, 0, 0, 'p', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'MarkerSize', 16);
    
    title(sprintf('\\textbf{Funzione di Lyapunov $V(x_k)$ (Partenza Reale)}'), 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('Deviazione $\Delta T_1$ [K]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Deviazione $\Delta T_2$ [K]', 'Interpreter', 'latex', 'FontSize', 12);
    zlabel('Energia $V(x_k) = x_k^T P_d x_k$', 'Interpreter', 'latex', 'FontSize', 12);
    
    legend([h_surf, h_line, h_start3, h_end], ...
           {'Superficie $V(x_k)$', 'Evoluzione LQR', 'Partenza', 'Origine'}, ...
           'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
    
    zlim([0, max(max_v_traj * 1.2, max(V_contour(:)))]); 
    xlim([-t_lim, t_lim]);
    ylim([-t_lim, t_lim]);
    view(-35, 30);
end