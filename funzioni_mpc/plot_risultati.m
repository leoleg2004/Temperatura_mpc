function plot_risultati(t_sim, storia_x, storia_u, U_min, U_max, x_ref, u_ref, Ts)
    if nargin < 8
        Ts = 120;
    end
    if nargin < 6
        x_ref = zeros(6,1);
        u_ref = zeros(3,1);
    end

    % PLOT_RISULTATI Disegna i grafici dell'evoluzione di stati e attuatori
    % Stati: T1, T2, T3 (Temperature) e Q1, Q2, Q3 (Potenze dei termosifoni)
    
    time_minutes = (0:t_sim) * Ts / 60;
    time_minutes_u = (0:t_sim-1) * Ts / 60;
    
    % Creazione di un'unica finestra globale ordinata
    figure('Name', 'Risultati MPC: Panoramica Globale', 'Color', 'w', 'Position', [100 100 1400 800]);
    
    % ================= COLONNA 1: TEMPERATURE =================
    subplot(3,3,1); 
    plot(time_minutes, storia_x(1,:), '-b', 'LineWidth', 1.5); hold on; 
    yline(x_ref(1), 'r--', 'LineWidth', 1.2); 
    title('Temperatura Stanza 1 (T_1) [K]'); grid on; ylabel('T [K]');
    ylim([x_ref(1)-6, x_ref(1)+3]);
    legend({'Traiettoria', 'Target'}, 'Location', 'best');
    
    subplot(3,3,4); 
    plot(time_minutes, storia_x(2,:), '-m', 'LineWidth', 1.5); hold on; 
    yline(x_ref(2), 'r--', 'LineWidth', 1.2); 
    title('Temperatura Stanza 2 (T_2) [K]'); grid on; ylabel('T [K]');
    ylim([x_ref(2)-6, x_ref(2)+3]);
    
    subplot(3,3,7); 
    plot(time_minutes, storia_x(3,:), '-c', 'LineWidth', 1.5); hold on; 
    yline(x_ref(3), 'r--', 'LineWidth', 1.2); 
    title('Temperatura Stanza 3 (T_3) [K]'); grid on; ylabel('T [K]'); xlabel('Tempo [min]');
    ylim([x_ref(3)-6, x_ref(3)+3]);

    % ================= COLONNA 2: CALORI EROGATI (STATI Q) =================
    subplot(3,3,2); 
    plot(time_minutes, storia_x(4,:), '-b', 'LineWidth', 1.5); hold on; 
    yline(x_ref(4), 'r--', 'LineWidth', 1.2); 
    title('Calore Termosifone 1 (Q_1) [W]'); grid on; ylabel('Q [W]');
    ylim([U_min(1)-10, U_max(1)+10]);
    legend({'Traiettoria', 'Target'}, 'Location', 'best');
    
    subplot(3,3,5); 
    plot(time_minutes, storia_x(5,:), '-m', 'LineWidth', 1.5); hold on; 
    yline(x_ref(5), 'r--', 'LineWidth', 1.2); 
    title('Calore Termosifone 2 (Q_2) [W]'); grid on; ylabel('Q [W]');
    ylim([U_min(2)-10, U_max(2)+10]);
    
    subplot(3,3,8); 
    plot(time_minutes, storia_x(6,:), '-c', 'LineWidth', 1.5); hold on; 
    yline(x_ref(6), 'r--', 'LineWidth', 1.2); 
    title('Calore Termosifone 3 (Q_3) [W]'); grid on; ylabel('Q [W]'); xlabel('Tempo [min]');
    ylim([U_min(3)-10, U_max(3)+10]);

    % ================= COLONNA 3: SFORZO ATTUATORI (RIFERIMENTI Q_r) =================
    subplot(3,3,3); stairs(time_minutes_u, storia_u(1,:), '-g', 'LineWidth', 1.5); hold on;
    yline(U_max(1), 'k--'); yline(U_min(1), 'k--'); yline(u_ref(1), 'r--', 'LineWidth', 1.2); 
    title(' Comando Ottimale Q_{1,r} [W]'); grid on; ylim([U_min(1)-10, U_max(1)+10]); ylabel('Potenza [W]');
    legend({'Comando MPC', 'Limiti', '', 'Target'}, 'Location', 'best');
    
    subplot(3,3,6); stairs(time_minutes_u, storia_u(2,:), '-g', 'LineWidth', 1.5); hold on;
    yline(U_max(2), 'k--'); yline(U_min(2), 'k--'); yline(u_ref(2), 'r--', 'LineWidth', 1.2);
    title(' Comando Ottimale Q_{2,r} [W]'); grid on; ylim([U_min(2)-10, U_max(2)+10]); ylabel('Potenza [W]');
    
    subplot(3,3,9); stairs(time_minutes_u, storia_u(3,:), '-g', 'LineWidth', 1.5); hold on;
    yline(U_max(3), 'k--'); yline(U_min(3), 'k--'); yline(u_ref(3), 'r--', 'LineWidth', 1.2);
    title(' Comando Ottimale Q_{3,r} [W]'); grid on; ylim([U_min(3)-10, U_max(3)+10]); ylabel('Potenza [W]'); xlabel('Tempo [min]');
end
