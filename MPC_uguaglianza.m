% =========================================================================
% Tesi Triennale - MPC Appartamento 3 Stanze
% FORMULAZIONE 2: VINCOLO TERMINALE DI UGUAGLIANZA
% =========================================================================

% Inizializzazione Path
[project_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(project_dir));
cd(project_dir);
disp('Path e cartella di lavoro inizializzati in automatico!');

%% 1. Definizione del Sistema (TEMPO CONTINUO)
nx = 6; % Stati: [T1, T2, T3, Q1, Q2, Q3]'
nu = 3; % Ingressi: [Q1_r, Q2_r, Q3_r]'
Ts = 120; % Tempo di campionamento in secondi

x_ref = [289; 289; 289; 100; 100; 100];
u_ref = [100; 100; 100];

disp('Ricerca dell''Operating Point ed estrazione delle matrici da Simulink...');
mdl = 'Stanza';
load_system(mdl);

opspec = operspec(mdl);
if ~isempty(opspec.States)
    opspec.States(1).x = x_ref;
    opspec.States(1).Known = [1; 1; 1; 0; 0; 0];
end
if ~isempty(opspec.Inputs)
    opspec.Inputs(1).u = u_ref;
    opspec.Inputs(1).Known = zeros(nu, 1);
end

opt = findopOptions('DisplayReport', 'off');
op = findop(mdl, opspec, opt);

x_ref = op.States(1).x;
if ~isempty(op.Inputs)
    u_ref = op.Inputs(1).u;
end

sys_c = linearize(mdl, op);
Ac = sys_c.A;
Bc = sys_c.B;

disp('Conversione del modello da Continuo a Discreto...');
[Ad, Bd] = discretizza_modello(Ac, Bc, nx, nu, Ts);

% -------------------------------------------------------------------------
% VERIFICA PROPRIETA' STRUTTURALI (Raggiungibilità e Osservabilità)
% -------------------------------------------------------------------------
run('Analisi_di_Sistema/verifica_raggiungibilita.m');

%% 2. Progetto Pesi Q e R
disp('Progetto LQR per estrarre i pesi Q e R (e matrice A_cl per Lyapunov)...');
[~, P, Q, R, A_cl] = progetta_LQR_discreto(Ad, Bd);

%% 3. Vincoli Fisici (Ampiezza e Rateo)
disp('Impostazione dei vincoli fisici (Ampiezza e Rateo)...');
U_min = [0; 0; 0];  
U_max = [150; 150; 150]; 


X_min = [282.5; 282.5; 282.5; 0; 0; 0];
X_max = [350; 350; 350; 150; 150; 150]; 

[Fx, fx, Fu, fu, dU_min, dU_max, Gx, gx] = imposta_vincoli(X_min, X_max, U_min, U_max, Ts);

%% 5. Setup Problema MPC con VINCOLO DI UGUAGLIANZA
% Il vincolo di uguaglianza costringe il sistema ad arrivare al 
% target in ESATTAMENTE N passi. Se N è piccolo (es. 6), per un sistema termico 
% lento è FISICAMENTE IMPOSSIBILE e quadprog restituirà INFEASIBLE.
% Per questo motivo, con l'uguaglianza è quasi obbligatorio usare un N molto grande!
N = 15; %devo mettere olti piu passi per aggiungere il punot di equilibrio
disp(['Setup MPC con Vincolo Terminale di Uguaglianza (x_N = x_ref) con N = ', num2str(N)]);

mpc_prob = setup_mpc_uguaglianza(N, nx, nu, Ad, Bd, Q, R, U_min, U_max, Gx, gx, x_ref, u_ref);

%% 6. Simulazione MPC Completa 
disp('--- Avvio Ottimizzazione e Simulazione MPC ---');
x_iniziale = [284; 285; 284; 0; 10; 0];   
t_sim = 150; 

[storia_x_eq, storia_u_eq, storia_costo_eq] = simula_mpc(mpc_prob, x_iniziale, t_sim, Ad, Bd, dU_max, x_ref, u_ref, '(Uguaglianza)');
disp('Ottimizzazione Riuscita!');

%% 7. Grafici 
plot_risultati(t_sim, storia_x_eq, storia_u_eq, U_min, U_max, x_ref, u_ref, Ts);

% Plot del Ritratto di Fase e della Funzione di Lyapunov (Partenza Reale)
disp('Generazione plot del Ritratto di Fase e Funzione di Lyapunov...');
plot_lyapunov_discrete(P, A_cl, Ts, x_iniziale, x_ref);

%% 9. Plot Funzionale di Costo 3D
disp('Generazione plot del Funzionale di Costo 3D...');
plot_mpc_cost_3d(mpc_prob, Ad, dU_max, storia_x_eq, storia_costo_eq, Ts);

disp('NOTA SUL CONTROLLABLE SET:');
disp('Con il vincolo terminale di uguaglianza, il set invariante O_inf è un SINGOLO PUNTO (x_ref).');
disp('Il Controllable set N-Step è estremamente difficile da plottare con N=100 e si omette in questa formulazione.');

%% 11. Calcolo dei passi necessari per la convergenza
tolleranza = 0.05; % Tolleranza di 0.05 Kelvin rispetto al target
passo_convergenza = -1;

for t = 1:size(storia_x_eq, 2)
    % Controlliamo se da questo passo in poi l'errore di TUTTE le stanze resta sotto la tolleranza
    errore_futuro_max = max(max(abs(storia_x_eq(1:3, t:end) - x_ref(1:3))));
    
    if errore_futuro_max <= tolleranza
        passo_convergenza = t - 1; % -1 perché t=1 rappresenta l'istante 0
        break;
    end
end

fprintf('\n======================================================\n');
if passo_convergenza >= 0
    fprintf('CONVERGENZA RAGGIUNTA (Tolleranza %.2f K)\n', tolleranza);
    fprintf('Passi necessari dall''MPC Uguaglianza: %d passi\n', passo_convergenza);
    fprintf('Tempo fisico di assestamento: %.1f minuti\n', (passo_convergenza * Ts)/60);
else
    disp('Convergenza NON raggiunta entro la fine della simulazione!');
    disp('Prova ad aumentare t_sim per dare più tempo al sistema.');
end
fprintf('======================================================\n');
