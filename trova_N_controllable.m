[project_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(project_dir));
cd(project_dir);
nx = 6; nu = 3; Ts = 120;
x_ref = [289; 289; 289; 100; 100; 100];
u_ref = [100; 100; 100];
mdl = 'Stanza';
load_system(mdl);
opspec = operspec(mdl);
opspec.States(1).x = x_ref; opspec.States(1).Known = [1; 1; 1; 0; 0; 0];
opspec.Inputs(1).u = u_ref; opspec.Inputs(1).Known = zeros(nu, 1);
opt = findopOptions('DisplayReport', 'off');
op = findop(mdl, opspec, opt);
sys_c = linearize(mdl, op);
Ac = sys_c.A; Bc = sys_c.B;
[Ad, Bd] = discretizza_modello(Ac, Bc, nx, nu, Ts);
[~, P, Q, R, A_cl] = progetta_LQR_discreto(Ad, Bd);
U_min = [0; 0; 0]; U_max = [150; 150; 150];
X_min = [282.5; 282.5; 282.5; 0; 0; 0]; X_max = [320; 320; 320; 150; 150; 150];
[Fx, fx, Fu, fu, dU_min, dU_max, Gx, gx] = imposta_vincoli(X_min, X_max, U_min, U_max, Ts);
[G_inf, g_inf] = cis(Ad, Bd, x_ref, u_ref, Fx, fx, Fu, fu, Q, R);

disp('Definizione del sistema MPT3...');
sys = LTISystem('A', Ad, 'B', Bd);
sys.x.min = X_min - x_ref;
sys.x.max = X_max - x_ref;
sys.u.min = U_min - u_ref;
sys.u.max = U_max - u_ref;
target = Polyhedron(G_inf, g_inf - G_inf*x_ref);

disp('Calcolo iterativo del Controllable Set all''indietro...');
S = target;
max_N = 100;
converged_N = -1;
for k = 1:max_N
    S_prev = S;
    S = sys.reachableSet('X', S_prev, 'direction', 'backward', 'N', 1);
    S.minHRep();
    if S == S_prev
        converged_N = k - 1;
        break;
    end
    if mod(k, 5) == 0
        fprintf('Iterazione N = %d completata...\n', k);
    end
end

if converged_N >= 0
    fprintf('\n>>> CONVERGENZA TROVATA! Il Controllable Set smette di crescere a N = %d passi! <<<\n', converged_N);
else
    fprintf('\nNessuna convergenza dopo %d passi.\n', max_N);
end
