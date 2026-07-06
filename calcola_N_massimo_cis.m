% Script per calcolare a quale N il Controllable Set smette di crescere
disp('Avvio algoritmo di ricerca per N massimo del Controllable Set...');

% Inizializzazione
startup_project();
nx = 6; nu = 3; Ts = 120;
x_ref = [289; 289; 289; 100; 100; 100];
u_ref = [100; 100; 100];

% Carica modello e linearizza
mdl = 'Stanza'; load_system(mdl);
opspec = operspec(mdl);
opspec.States(1).x = x_ref; opspec.States(1).Known = [1; 1; 1; 0; 0; 0];
opspec.Inputs(1).u = u_ref; opspec.Inputs(1).Known = zeros(nu, 1);
opt = findopOptions('DisplayReport', 'off'); op = findop(mdl, opspec, opt);
sys_c = linearize(mdl, op);
Ac = sys_c.A; Bc = sys_c.B;
[Ad, Bd] = discretizza_modello(Ac, Bc, nx, nu, Ts);

% LQR e Vincoli
[~, P, Q, R, A_cl] = progetta_LQR_discreto(Ad, Bd);
U_min = [0; 0; 0]; U_max = [150; 150; 150];
X_min = [282.5; 282.5; 282.5; 0; 0; 0]; X_max = [320; 320; 320; 150; 150; 150];
[Fx, fx, Fu, fu, dU_min, dU_max, Gx, gx] = imposta_vincoli(X_min, X_max, U_min, U_max, Ts);

% Trova O_inf
[G_inf, g_inf] = cis(Ad, Bd, x_ref, u_ref, Fx, fx, Fu, fu, Q, R);
g_inf_zero = g_inf - G_inf*x_ref;
Hx_mpt = Fx; hx_mpt = fx - Fx*x_ref;
Hu_mpt = [eye(nu); -eye(nu)]; hu_mpt = [U_max - u_ref; u_ref - U_min];

% Calcolo Iterativo
disp('Calcolo espansione all''indietro... attendere...');
N_max = 100;
P_prev = [];
convergenza_N = -1;

for N_test = 1:N_max
    % Calcola la slice delle temperature (T1, T2, T3) bypassando il 6D
    P_curr = controllable_set_slice(Hx_mpt, hx_mpt, Hu_mpt, hu_mpt, G_inf, g_inf_zero, Ad, Bd, N_test, [4,5,6], [0;0;0]);
    P_curr.minVRep(); % Forza il calcolo dei vertici per pulizia geometrica
    
    if ~isempty(P_prev)
        % Se il poliedro attuale è uguale al precedente, abbiamo trovato il limite!
        if P_curr == P_prev
            convergenza_N = N_test - 1;
            break;
        end
    end
    P_prev = P_curr;
    if mod(N_test, 5) == 0
        fprintf('Passo N = %d calcolato...\n', N_test);
    end
end

fprintf('\n*******************************************************\n');
if convergenza_N > 0
    fprintf('BINGO! Il Controllable Set ha smesso di crescere a N = %d!\n', convergenza_N);
    fprintf('A N = %d, il poliedro ha "sbattuto" contro i vincoli termici.\n', convergenza_N);
else
    fprintf('Il poliedro continua a crescere fino a N = %d.\n', N_max);
end
fprintf('*******************************************************\n');
