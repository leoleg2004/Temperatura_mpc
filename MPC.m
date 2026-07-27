% =========================================================================
% Tesi Triennale - MPC Longitudinale F-16
% Simulazione con Conversione Continuo-Discreto (c2d) e Radianti
% Ing. Leggeri Leonardo
% =========================================================================

% Inizializzazione Path
startup_project();

%% 1. Definizione del Sistema (TEMPO CONTINUO)
nx = 6; % Stati: [T1, T2, T3, Q1, Q2, Q3]'
nu = 3; % Ingressi: [Q1_r, Q2_r, Q3_r]'
Ts = 120; % Tempo di campionamento in secondi (aumentato per dinamica termica lenta)

% Riferimento (Equilibrio desiderato dal progetto)
x_ref = [289; 289; 289; 100; 100; 100];
u_ref = [100; 100; 100];

% =========================================================================
% Estrazione del modello linearizzato esatto tramite Simulink
% =========================================================================
disp('Ricerca dell''Operating Point ed estrazione delle matrici da Simulink...');
mdl = 'Stanza';
load_system(mdl);

% Specifica del punto di funzionamento
opspec = operspec(mdl);

% Per trovare l'equilibrio termico perfetto, fissiamo le Temperature (289 K)
% ma lasciamo LIBERI gli stati dei termosifoni (Q1, Q2, Q3) e gli ingressi (u)
if ~isempty(opspec.States)
    opspec.States(1).x = x_ref;
    % [1;1;1] blocca le T. [0;0;0] lascia liberi i calori Q.
    opspec.States(1).Known = [1; 1; 1; 0; 0; 0];
end

if ~isempty(opspec.Inputs)
    opspec.Inputs(1).u = u_ref;
    opspec.Inputs(1).Known = zeros(nu, 1); % Ingressi completamente liberi
end

% Trova il vero punto operativo di equilibrio
opt = findopOptions('DisplayReport', 'off');
op = findop(mdl, opspec, opt);

% Aggiorniamo i nostri target x_ref e u_ref con i valori esatti di 
% equilibrio termodinamico trovati da Simulink (es. ~99W invece di 100W)
x_ref = op.States(1).x;
if ~isempty(op.Inputs)
    u_ref = op.Inputs(1).u;
end

% Estrai le matrici Ac e Bc
sys_c = linearize(mdl, op);
Ac = sys_c.A;
Bc = sys_c.B;

if isempty(Bc)
    warning('Attenzione: la matrice Bc estratta da Simulink è vuota! Assicurati di aver collegato correttamente una Inport al blocco Matlab Function.');
end

% Manteniamo l'offset affine calcolato analiticamente
f_eq = dinamica_casa(x_ref, u_ref);

disp('Conversione del modello da Continuo a Discreto...');
[Ad, Bd] = discretizza_modello(Ac, Bc, nx, nu, Ts);

% -------------------------------------------------------------------------
% VERIFICA PROPRIETA' STRUTTURALI (Raggiungibilità e Osservabilità)
% -------------------------------------------------------------------------
run('Analisi_di_Sistema/verifica_raggiungibilita.m');

%% 2. Progetto LQR tramite funzione dedicata
% =========================================================================
% NOTA SUI PESI Q e R:
% Se vuoi modificare l'aggressività del controllore (pesi Q e R),
% apri il file "funzioni_lqr/progetta_LQR_discreto.m" e modificali lì dentro!
% =========================================================================
disp('Progetto LQR discreto tramite funzione dedicata...');
[K, P, Q, R, A_cl] = progetta_LQR_discreto(Ad, Bd);

%% 3. Vincoli Fisici (Ampiezza e Rateo)
disp('Impostazione dei vincoli fisici (Ampiezza e Rateo)...');
% Ingressi: [Q1_r, Q2_r, Q3_r] (Potenze in Watt)
U_min = [0; 0; 0];  
U_max = [150; 150; 150]; 

% Stati: T1(K), T2(K), T3(K), Q1(W), Q2(W), Q3(W)
X_min = [282.5; 282.5; 282.5; 0; 0; 0];
X_max = [320; 320; 320; 150; 150; 150]; 

[Fx, fx, Fu, fu, dU_min, dU_max, Gx, gx] = imposta_vincoli(X_min, X_max, U_min, U_max, Ts);

% Nota: dato che il modello è stato linearizzato esattamente attorno a x_ref, u_ref,
% esso costituirà un ottimo punto di equilibrio locale.

%% 4. Calcolo Control Invariant Set e Plot
disp('--- CALCOLO del Control Invariant Set ---');
[G_inf, g_inf] = cis(Ad, Bd, x_ref, u_ref, Fx, fx, Fu, fu, Q, R);

%% 5. Setup Problema MPC 
N = 5; % Orizzonte predittivo 
mpc_prob = setup_mpc(N, nx, nu, Ad, Bd, Q, P, R, U_min, U_max, Gx, gx, G_inf, g_inf, x_ref, u_ref);

%% 6. Simulazione MPC Completa 
disp('--- Avvio Ottimizzazione e Simulazione MPC ---');
% Ordine: [T1; T2; T3; Q1; Q2; Q3]
x_iniziale = [284;  % T1
              285;  % T2
              284;  % T3
              0;    % Q1
              10;   % Q2
              0];   % Q3
t_sim = 120; 

[storia_x, storia_u, storia_costo] = simula_mpc(mpc_prob, x_iniziale, t_sim, Ad, Bd, dU_max, x_ref, u_ref);
disp('Ottimizzazione Riuscita!');

%% 7. Grafici 
% Disegna il Control Invariant Set in 3D con il punto di arrivo reale
plot_cis(G_inf, g_inf, x_ref, storia_x);

% Plot nel dominio del tempo
plot_risultati(t_sim, storia_x, storia_u, U_min, U_max, x_ref, u_ref, Ts);

% Plot del Ritratto di Fase e della Funzione di Lyapunov (Partenza Reale)
disp('Generazione plot del Ritratto di Fase e Funzione di Lyapunov...');
plot_lyapunov_discrete(P, A_cl, Ts, x_iniziale, x_ref);

%% 9. Plot Funzionale di Costo 3D
disp('Generazione plot del Funzionale di Costo 3D...');
plot_mpc_cost_3d(mpc_prob, Ad, dU_max, storia_x, storia_costo, Ts);

%% 10. Calcolo e Plot N-Step Controllable Set (tramite MPT3)
% Controlliamo se la classe Polyhedron (MPT3) esiste
if exist('Polyhedron', 'class') == 8
    disp('Calcolo N-Step Controllable Set (richiede MPT3)...');
    N_steps_ctrl = N; % Usa automaticamente l'orizzonte predittivo dell'MPC
    disp(['Calcolo Controllable Set per N = ', num2str(N_steps_ctrl), ' passi. Attenzione: MPT3 potrebbe impiegare molto tempo!']);
    
    Hx_mpt = Fx;
    hx_mpt = fx - Fx*x_ref;
    
    % Delta Ingressi ammissibili
    Hu_mpt = [eye(nu); -eye(nu)]; 
    hu_mpt = [U_max - u_ref; u_ref - U_min];
    
    % Il target è il Control Invariant Set che abbiamo già calcolato prima
    % Trasliamo g_inf allo zero
    g_inf_zero = g_inf - G_inf*x_ref;
    
    % Invece di calcolare il poliedro 6D completo (che bloccherebbe il PC per N alto),
    % calcoliamo DIRETTAMENTE le fette 3D richieste bypassando la proiezione 6D!
    disp('Plot del N-Step Controllable Set in corso (Algoritmo Super-Veloce 3D)...');
    
    % Dato che lo stato è a 6 dimensioni [T1, T2, T3, Q1, Q2, Q3],
    % "affettiamo" il poliedro bloccando gli stati dei termosifoni (Q) a 0.
    poly_slice_T = controllable_set_slice(Hx_mpt, hx_mpt, Hu_mpt, hu_mpt, G_inf, g_inf_zero, Ad, Bd, N_steps_ctrl, [4,5,6], [0;0;0]); 
    poly_slice_T = poly_slice_T + x_ref(1:3); % Trasliamo sulle vere temperature target
    
    % Calcoliamo una trasparenza dinamica: più passi fai, più è trasparente!
    alpha_val = max(0.1, 0.5 - (N_steps_ctrl - 1) * 0.1); 
    
    fig_cis_temp = findobj('type', 'figure', 'name', 'CIS 3D: Temperature (T1, T2, T3)');
    if ~isempty(fig_cis_temp)
        figure(fig_cis_temp(1)); hold on;
        % Disegniamo il Controllable Set (Giallo) con trasparenza dinamica
        h_ctrl_plot_T = poly_slice_T.plot('Alpha', alpha_val, 'Color', 'y');
        
        lgd = legend;
        if ~isempty(lgd)
            lgd.String{end} = sprintf('N-Step Controllable Set (N=%d)', N_steps_ctrl);
        end
        title(sprintf('\\textbf{CIS $\\mathcal{O}_\\infty$ e %d-Step Controllable Set}', N_steps_ctrl), 'Interpreter', 'latex', 'FontSize', 14);
    end
    
    % =========================================================
    % 2. Plot Controllable Set per i CALORI
    % =========================================================
    % "affettiamo" bloccando gli stati delle temperature (T) a 0.
    poly_slice_Q = controllable_set_slice(Hx_mpt, hx_mpt, Hu_mpt, hu_mpt, G_inf, g_inf_zero, Ad, Bd, N_steps_ctrl, [1,2,3], [0;0;0]); 
    poly_slice_Q = poly_slice_Q + x_ref(4:6); % Trasliamo sui veri calori target
    
    fig_cis_calori = findobj('type', 'figure', 'name', 'CIS 3D: Calori (Q1, Q2, Q3)');
    if ~isempty(fig_cis_calori)
        figure(fig_cis_calori(1)); hold on;
        % Disegniamo il Controllable Set (Giallo) con trasparenza dinamica
        h_ctrl_plot_Q = poly_slice_Q.plot('Alpha', alpha_val, 'Color', 'y');
        
        lgd2 = legend;
        if ~isempty(lgd2)
            lgd2.String{end} = sprintf('N-Step Controllable Set (N=%d)', N_steps_ctrl);
        end
        title(sprintf('\\textbf{CIS $\\mathcal{O}_\\infty$ e %d-Step Controllable Set}', N_steps_ctrl), 'Interpreter', 'latex', 'FontSize', 14);
    end
else
    disp('ATTENZIONE: Il toolbox MPT3 non è installato in questo MATLAB. Impossibile calcolare il Controllable Set in N-Step. Per visualizzarlo, scarica e installa MPT3 (https://www.mpt3.org).');
end
%% 11. Calcolo dei passi necessari per la convergenza
tolleranza = 0.05; % Tolleranza di 0.05 Kelvin rispetto al target
passo_convergenza = -1;

for t = 1:size(storia_x, 2)
    % Controlliamo se da questo passo in poi l'errore di TUTTE le stanze resta sotto la tolleranza
    errore_futuro_max = max(max(abs(storia_x(1:3, t:end) - x_ref(1:3))));
    
    if errore_futuro_max <= tolleranza
        passo_convergenza = t - 1; % -1 perché t=1 rappresenta l'istante 0
        break;
    end
end

fprintf('\n======================================================\n');
if passo_convergenza >= 0
    fprintf('CONVERGENZA RAGGIUNTA (Tolleranza %.2f K)\n', tolleranza);
    fprintf('Passi necessari dall''MPC: %d passi\n', passo_convergenza);
    fprintf('Tempo fisico di assestamento: %.1f minuti\n', (passo_convergenza * Ts)/60);
else
    disp('Convergenza NON raggiunta entro la fine della simulazione!');
    disp('Prova ad aumentare t_sim per dare più tempo al sistema.');
end
fprintf('======================================================\n');
