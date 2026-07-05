function [K, P, Q, R, A_cl] = progetta_LQR_discreto(A, B)
    % =====================================================================
    % Calcolo pesi LQR coerenti con il modello Casa (3 stanze)
    % Ordine stati: [T1; T2; T3; Q1; Q2; Q3]
    % Ordine input: [Q1_r; Q2_r; Q3_r]
    % =====================================================================
    
    % Massimi scostamenti tollerabili (Regola di Bryson)
    max_T = 2;    % 2 Kelvin di errore max tollerato sulle temperature
    max_Q = 20;   % 20 Watt di scostamento sui calori erogati (stati attuatore)
    
    Q = diag([1/max_T^2, 1/max_T^2, 1/max_T^2, 1/max_Q^2, 1/max_Q^2, 1/max_Q^2]);
    
    % Massimi scostamenti tollerabili per gli ingressi (riferimenti di potenza)
    max_Qr = 50;  % 50 Watt di variazione sui riferimenti degli attuatori
    
    R = diag([1/max_Qr^2, 1/max_Qr^2, 1/max_Qr^2]);
    
    % Fattori di tuning
    rho_q = 1; 
    rho_r = 100; 
    
    Q = rho_q * Q;
    R = rho_r * R;
    
    % Sintesi LQR
    [K, P, ~] = dlqr(A, B, Q, R);
    A_cl = A - B*K;
end
