% =========================================================================
% ANALISI DI RAGGIUNGIBILITA' E OSSERVABILITA' 
% Progetto MPC - Appartamento 3 Stanze
% =========================================================================

% Assicuriamoci che le matrici del modello siano caricate in memoria.
% Se non lo sono, eseguiamo lo script principale fino al punto di linearizzazione.
if ~exist('Ad', 'var') || ~exist('Bd', 'var')
    disp('Matrici Ad e Bd non trovate nel workspace. Eseguo l''inizializzazione del modello...');
    % Non eseguiamo tutto MPC.m per non lanciare ottimizzazioni lunghe,
    % ma avviamo la parte di setup.
    disp('Avvia prima MPC.m per avere le matrici di stato a disposizione!');
    return;
end

%% Calcolo della Raggiungibilità / Controllabilità
n = size(Ad, 1);       

% 1. Calcolo della matrice di Raggiungibilità (Controllabilità in discreto)
R = ctrb(Ad, Bd);

% 2. Calcolo del rango della matrice
rango_R = rank(R);    

disp('--- Analisi di Raggiungibilità/Controllabilità ---');
disp('Matrice di Raggiungibilità R (Prime 6 colonne):'); 
disp(R(:, 1:6)); % Stampiamo solo le prime colonne per non inondare lo schermo

% 3. Verifica della completa raggiungibilità
if rango_R == n
    disp(['✅ Il sistema termico è COMPLETAMENTE raggiungibile (Rango R = ', num2str(rango_R), ').']);
else
    disp('❌ Il sistema termico NON è completamente raggiungibile (rango incompleto).');
    disp(['Il rango della matrice R è ', num2str(rango_R), ' invece di ', num2str(n), '.']);
end

%% Calcolo dell'Osservabilità (ipotizzando tutti gli stati misurabili)
% Nel nostro progetto controlliamo interamente lo stato, quindi supponiamo C = eye(nx)
Cd = eye(n);

% 1. Calcolo della matrice di Osservabilità
O = obsv(Ad, Cd);

% 2. Calcolo del rango della matrice
rango_O = rank(O);    

disp(' '); 
disp('--- Analisi di Osservabilità ---');
disp('Matrice di Osservabilità O (Prime 6 righe):'); 
disp(O(1:6, :)); 

% 3. Verifica della completa osservabilità
if rango_O == n
    disp(['✅ Il sistema termico è COMPLETAMENTE osservabile (Rango O = ', num2str(rango_O), ').']);
else
    disp('❌ Il sistema termico NON è completamente osservabile (rango incompleto).');
    disp(['Il rango della matrice O è ', num2str(rango_O), ' invece di ', num2str(n), '.']);
end
