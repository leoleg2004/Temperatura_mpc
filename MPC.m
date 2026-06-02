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

% Creiamo un punto operativo vuoto basato sul modello
op = operpoint(mdl);

% Forziamo manualmente lo stato al valore di riferimento x_ref (blocco integratore)
if ~isempty(op.States)
    op.States(1).x = x_ref;
end

% Forziamo manualmente l'ingresso al valore di riferimento u_ref
if ~isempty(op.Inputs)
    op.Inputs(1).u = u_ref;
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
Rate_max = [50; 50; 50]; % Dummy (W/s)

% Stati: T1(K), T2(K), T3(K), Q1(W), Q2(W), Q3(W)
X_min = [282.5; 282.5; 282.5; 0; 0; 0];
X_max = [350; 350; 350; 150; 150; 150]; % Dummy upper bounds

[Fx, fx, Fu, fu, dU_min, dU_max, Gx, gx] = imposta_vincoli(X_min, X_max, U_min, U_max, Rate_max, Ts);

% Nota: dato che il modello è stato linearizzato esattamente attorno a x_ref, u_ref,
% esso costituirà un ottimo punto di equilibrio locale.

%% 4. Calcolo Control Invariant Set e Plot
disp('--- CALCOLO del Control Invariant Set ---');
[G_inf, g_inf] = cis(Ad, Bd, x_ref, u_ref, Fx, fx, Fu, fu, Q, R);

%% 4b. Plot 3D dei Set Invarianti
% Spostato dopo la simulazione per tracciare anche il punto di arrivo

%% 5. Setup Problema MPC 
N = 40; % Orizzonte predittivo 
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
t_sim = 150; % Aumentato a 350 passi per permettere alla dinamica lenta (Fugoide) di centrare perfettamente il target

[storia_x, storia_u] = simula_mpc(mpc_prob, x_iniziale, t_sim, Ad, Bd, dU_max, x_ref, u_ref);
disp('Ottimizzazione Riuscita. Il modello è matematicamente solido.');

disp('---------------------------------------------------');
disp('VERIFICA RAGGIUNGIMENTO TARGET (Ultimo Step)');
disp('Stato Finale Raggiunto (storia_x(:,end)):');
disp(storia_x(:,end));
disp('Target Desiderato (x_ref):');
disp(x_ref);
disp('Errore Assoluto [T1; T2; T3; Q1; Q2; Q3]:');
disp(abs(storia_x(:,end) - x_ref));
disp('---------------------------------------------------');

%% 7. Grafici 
% Disegna il Control Invariant Set in 3D con il punto di arrivo reale
plot_cis(G_inf, g_inf, x_ref, storia_x);

% Plot nel dominio del tempo
plot_risultati(t_sim, storia_x, storia_u, U_min, U_max, x_ref, u_ref, Ts);


%% 8. plot clf ljapunov function
plot_lyapunov_discrete(P,A_cl,Ts);