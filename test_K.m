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
[K, P, Q, R, A_cl] = progetta_LQR_discreto(Ad, Bd);

x_iniziale = [288.5; 288.5; 288.5; 100; 100; 100];
u_unconstrained = -K * (x_iniziale - x_ref) + u_ref;
disp('u_unconstrained = ');
disp(u_unconstrained);
