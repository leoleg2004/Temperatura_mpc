function plot_lyapunov_discrete(P_ds, A_cl_ds, Ts)
    % =========================================================================
    % PLOT_LYAPUNOV_DISCRETE - Funzione di Lyapunov 3D per sistemi a tempo discreto
    % Grafico migliorato: Traiettoria interpolata fluida + Campionamenti discreti
    % =========================================================================
    
    if nargin < 3
        Ts = 0.05; 
    end
    
    rad2deg = 180 / pi;
    
    % --- 1. DEFINIZIONE LIMITI VISIVI ---
    t_lim = 15;  % Kelvin (deviazione dal target)
    
    [T1_grid, T2_grid] = meshgrid(linspace(-t_lim, t_lim, 80), linspace(-t_lim, t_lim, 80));
    
    % --- 2. CALCOLO SUPERFICIE DI LYAPUNOV ---
    V_surf = zeros(size(T1_grid));
    for i = 1:size(T1_grid, 1)
        for j = 1:size(T1_grid, 2)
            stato_surf = [T1_grid(i,j); T2_grid(i,j); 0; 0; 0; 0];
            V_surf(i,j) = stato_surf' * P_ds * stato_surf;
        end
    end
    
    % --- 3. DINAMICA DISCRETA E CONDIZIONI INIZIALI ---
    N_steps = round(300 / Ts); % Simuliamo un po' di passi
    if N_steps < 10; N_steps = 10; end
    
    x0 = [ 5,  3, 0, 0, 0, 0;   
          -4, -6, 0, 0, 0, 0;   
           2, -5, 0, 0, 0, 0]';

    %% --- FIGURA 3D ---
    figure('Name', sprintf('Lyapunov 3D Discreto (Ts = %g s)', Ts), 'Color', 'w', 'Position', [150 150 850 650]);
    hold on; grid on;
    
    % 1. Disegno Superficie
    h_surf = surf(T1_grid, T2_grid, V_surf, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
    colormap jet;
    cb = colorbar;
    ylabel(cb, 'Energia $V(x_k)$', 'Interpreter', 'latex', 'FontSize', 12);
    
    colori = {'r', 'g', 'b'};
    max_v_traj = 0; 
    
    % 2. Simulazione e Disegno Traiettorie
    for i = 1:size(x0, 2)
        x_traj = zeros(6, N_steps);
        x_traj(:, 1) = x0(:, i);
        
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
        
        % --- MAGIA VISIVA: INTERPOLAZIONE FLUIDA DELLA TRAIETTORIA ---
        % Usiamo 'pchip' (interpolazione cubica che non crea finti rimbalzi) 
        % per generare 10 volte più punti e rendere la curva morbidissima
        t_discrete = 1:N_steps;
        t_fine = linspace(1, N_steps, N_steps * 10);
        
        T1_smooth = pchip(t_discrete, x_traj(1,:), t_fine);
        T2_smooth = pchip(t_discrete, x_traj(2,:), t_fine);
        V_smooth = pchip(t_discrete, V_traj, t_fine);
        
        % 1. Disegniamo prima la SCIA FLUIDA (linea continua)
        if i == 1
            h_line = plot3(T1_smooth, T2_smooth, V_smooth, '-', 'Color', colori{i}, 'LineWidth', 2);
        else
            plot3(T1_smooth, T2_smooth, V_smooth, '-', 'Color', colori{i}, 'LineWidth', 2);
        end
        
        % 2. Disegniamo sopra i PUNTI DISCRETI (solo pallini, senza linee)
        plot3(x_traj(1,:), x_traj(2,:), V_traj, 'o', ...
              'MarkerEdgeColor', colori{i}, 'MarkerFaceColor', 'w', 'MarkerSize', 4);
              
        % Partenza (Start)
        if i == 1
            h_start = plot3(x_traj(1,1), x_traj(2,1), V_traj(1), 'o', ...
                            'MarkerFaceColor', colori{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        else
            plot3(x_traj(1,1), x_traj(2,1), V_traj(1), 'o', ...
                  'MarkerFaceColor', colori{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        end
        
        % Arrivo (End)
        plot3(x_traj(1,end), x_traj(2,end), V_traj(end), 's', ...
              'MarkerFaceColor', colori{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    
    % --- ORIGINE ---
    h_end = plot3(0, 0, 0, 'p', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'MarkerSize', 16);
    
    % --- FORMATTAZIONE E LEGENDA ---
    title(sprintf('\\textbf{Funzione di Lyapunov $V(x_k)$ a Tempo Discreto ($T_s = %g$ s)}', Ts), 'Interpreter', 'latex', 'FontSize', 16);
    xlabel('Deviazione $\Delta T_1$ [K]', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Deviazione $\Delta T_2$ [K]', 'Interpreter', 'latex', 'FontSize', 12);
    zlabel('Energia $V(x_k) = x_k^T P_d x_k$', 'Interpreter', 'latex', 'FontSize', 12);
    
    legend([h_surf, h_line, h_start, h_end], ...
           {'Superficie $V(x_k)$', 'Evoluzione di Stato (Interpolata)', 'Partenza ($k=0$)', 'Origine ($k \to \infty$)'}, ...
           'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast');
    
    zlim([0, max_v_traj * 1.2]); 
    xlim([-t_lim, t_lim]);
    ylim([-t_lim, t_lim]);
    view(-35, 30);
end