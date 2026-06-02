function plot_cis(G_inf, g_inf, x_ref, storia_x)
    if nargin < 4
        storia_x = [];
    end
    if nargin < 3
        x_ref = zeros(6,1);
    end
    
    % PLOT_CIS Disegna in 3D le fette del Control Invariant Set per le 3 stanze
    
    %% PLOT 1: Fetta delle Temperature (Fissando i Calori ai valori di target)
    figure('Name', 'CIS 3D: Temperature (T1, T2, T3)', 'Color', 'w', 'Position', [50, 50, 800, 600]);
    hold on; grid on; view(3); 
    title(sprintf('Poliedro $\\mathcal{X}_f$ in 3D - Temperature\nFetta: $Q_1=%.1f, Q_2=%.1f, Q_3=%.1f$', x_ref(4), x_ref(5), x_ref(6)), 'Interpreter', 'latex', 'FontSize', 12);
    
    xlabel('T_1 [K]'); ylabel('T_2 [K]'); zlabel('T_3 [K]');
    
    % Griglia per T1, T2, T3 attorno al target
    [X_T1, X_T2, X_T3] = meshgrid(linspace(x_ref(1)-15, x_ref(1)+15, 30), ...
                                  linspace(x_ref(2)-15, x_ref(2)+15, 30), ...
                                  linspace(x_ref(3)-15, x_ref(3)+15, 30));
    X_T1_f = X_T1(:); X_T2_f = X_T2(:); X_T3_f = X_T3(:);
    
    % Fissiamo i calori (stati 4,5,6)
    X_Q1_f = x_ref(4) * ones(size(X_T1_f));
    X_Q2_f = x_ref(5) * ones(size(X_T1_f));
    X_Q3_f = x_ref(6) * ones(size(X_T1_f));
    
    Validi_inf_T = true(size(X_T1_f));
    for j = 1:size(G_inf, 1)
        Valore = G_inf(j,1)*X_T1_f + G_inf(j,2)*X_T2_f + G_inf(j,3)*X_T3_f + ...
                 G_inf(j,4)*X_Q1_f + G_inf(j,5)*X_Q2_f + G_inf(j,6)*X_Q3_f;
        Validi_inf_T = Validi_inf_T & (Valore <= g_inf(j));
    end
    
    PX_T1 = X_T1_f(Validi_inf_T); PY_T2 = X_T2_f(Validi_inf_T); PZ_T3 = X_T3_f(Validi_inf_T);
    
    if length(PX_T1) > 4
        try
            K_hull_T = convhull(PX_T1, PY_T2, PZ_T3);
            trisurf(K_hull_T, PX_T1, PY_T2, PZ_T3, 'FaceColor', 'r', 'FaceAlpha', 0.4, 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
        catch
            plot3(PX_T1, PY_T2, PZ_T3, 'r.', 'MarkerSize', 10);
        end
    else
        disp('ATTENZIONE: Nessun punto valido trovato per il plot di O_inf sulla fetta delle temperature.');
    end
    
    % Target
    plot3(x_ref(1), x_ref(2), x_ref(3), 'kX', 'MarkerSize', 12, 'LineWidth', 3);
    
    if ~isempty(storia_x)
        plot3(storia_x(1,:), storia_x(2,:), storia_x(3,:), '-m', 'LineWidth', 2);
        plot3(storia_x(1,:), storia_x(2,:), storia_x(3,:), 'mo', 'MarkerSize', 4, 'MarkerFaceColor', 'm');
        plot3(storia_x(1,1), storia_x(2,1), storia_x(3,1), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
        plot3(storia_x(1,end), storia_x(2,end), storia_x(3,end), 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g');
        legend({'CIS $\mathcal{O}_\infty$ (Temperature)', 'Target', 'Traiettoria MPC', 'Passi', 'Partenza', 'Arrivo'}, 'Location', 'best', 'Interpreter', 'latex');
    else
        legend({'CIS $\mathcal{O}_\infty$ (Temperature)', 'Target'}, 'Location', 'best', 'Interpreter', 'latex');
    end
    
    %% PLOT 2: Fetta dei Calori (Fissando le Temperature ai valori di target)
    figure('Name', 'CIS 3D: Calori (Q1, Q2, Q3)', 'Color', 'w', 'Position', [100, 100, 800, 600]);
    hold on; grid on; view(3); 
    title(sprintf('Poliedro $\\mathcal{X}_f$ in 3D - Calori\nFetta: $T_1=%.1f, T_2=%.1f, T_3=%.1f$', x_ref(1), x_ref(2), x_ref(3)), 'Interpreter', 'latex', 'FontSize', 12);
    
    xlabel('Q_1 [W]'); ylabel('Q_2 [W]'); zlabel('Q_3 [W]');
    
    % Griglia per Q1, Q2, Q3 attorno al target
    [X_Q1, X_Q2, X_Q3] = meshgrid(linspace(x_ref(4)-25, x_ref(4)+25, 30), ...
                                  linspace(x_ref(5)-25, x_ref(5)+25, 30), ...
                                  linspace(x_ref(6)-25, x_ref(6)+25, 30));
    X_Q1_f = X_Q1(:); X_Q2_f = X_Q2(:); X_Q3_f = X_Q3(:);
    
    % Fissiamo le temperature (stati 1,2,3)
    X_T1_f = x_ref(1) * ones(size(X_Q1_f));
    X_T2_f = x_ref(2) * ones(size(X_Q1_f));
    X_T3_f = x_ref(3) * ones(size(X_Q1_f));
    
    Validi_inf_Q = true(size(X_Q1_f));
    for j = 1:size(G_inf, 1)
        Valore = G_inf(j,1)*X_T1_f + G_inf(j,2)*X_T2_f + G_inf(j,3)*X_T3_f + ...
                 G_inf(j,4)*X_Q1_f + G_inf(j,5)*X_Q2_f + G_inf(j,6)*X_Q3_f;
        Validi_inf_Q = Validi_inf_Q & (Valore <= g_inf(j));
    end
    
    PX_Q1 = X_Q1_f(Validi_inf_Q); PY_Q2 = X_Q2_f(Validi_inf_Q); PZ_Q3 = X_Q3_f(Validi_inf_Q);
    
    if length(PX_Q1) > 4
        try
            K_hull_Q = convhull(PX_Q1, PY_Q2, PZ_Q3);
            trisurf(K_hull_Q, PX_Q1, PY_Q2, PZ_Q3, 'FaceColor', 'b', 'FaceAlpha', 0.4, 'EdgeColor', 'k', 'EdgeAlpha', 0.2);
        catch
            plot3(PX_Q1, PY_Q2, PZ_Q3, 'b.', 'MarkerSize', 10);
        end
    else
        disp('ATTENZIONE: Nessun punto valido trovato per il plot di O_inf sulla fetta dei calori.');
    end
    
    % Target
    plot3(x_ref(4), x_ref(5), x_ref(6), 'kX', 'MarkerSize', 12, 'LineWidth', 3);
    
    if ~isempty(storia_x)
        plot3(storia_x(4,:), storia_x(5,:), storia_x(6,:), '-m', 'LineWidth', 2);
        plot3(storia_x(4,:), storia_x(5,:), storia_x(6,:), 'mo', 'MarkerSize', 4, 'MarkerFaceColor', 'm');
        plot3(storia_x(4,1), storia_x(5,1), storia_x(6,1), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
        plot3(storia_x(4,end), storia_x(5,end), storia_x(6,end), 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g');
        legend({'CIS $\mathcal{O}_\infty$ (Calori)', 'Target', 'Traiettoria MPC', 'Passi', 'Partenza', 'Arrivo'}, 'Location', 'best', 'Interpreter', 'latex');
    else
        legend({'CIS $\mathcal{O}_\infty$ (Calori)', 'Target'}, 'Location', 'best', 'Interpreter', 'latex');
    end
end
