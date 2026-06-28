function [poly_slice_out] = controllable_set_slice(Hx, hx, Hu, hu, H_target, h_target, A, B, N, slice_dims, slice_vals)
% Calcola direttamente una SLICE 3D dell'N-Step Controllable Set, bypassando 
% il calcolo esponenziale in 6D. Questo riduce il tempo di calcolo da minuti/ore a decimi di secondo!

n = size(A, 1);
m = size(B, 2);

plot_dims = setdiff(1:n, slice_dims);
n_plot = length(plot_dims);

% Costruiamo le matrici estese per x_0 e U = [u_0; ...; u_{N-1}]
% Dimensione di U: N*m
A_x = [];
A_U = [];
b_tot = [];

% 1. Vincoli di stato: Hx * x_k <= hx per k = 0 ... N-1
for k = 0:N-1
    Phi_k = A^k;
    Gamma_k = zeros(n, N*m);
    for i = 0:k-1
        Gamma_k(:, i*m+1 : (i+1)*m) = A^(k-1-i) * B;
    end
    A_x = [A_x; Hx * Phi_k];
    A_U = [A_U; Hx * Gamma_k];
    b_tot = [b_tot; hx];
end

% 2. Vincoli di ingresso: Hu * u_k <= hu per k = 0 ... N-1
for k = 0:N-1
    Gamma_U = zeros(size(Hu, 1), N*m);
    Gamma_U(:, k*m+1 : (k+1)*m) = Hu;
    A_x = [A_x; zeros(size(Hu, 1), n)];
    A_U = [A_U; Gamma_U];
    b_tot = [b_tot; hu];
end

% 3. Vincolo Target: H_target * x_N <= h_target
Phi_N = A^N;
Gamma_N = zeros(n, N*m);
for i = 0:N-1
    Gamma_N(:, i*m+1 : (i+1)*m) = A^(N-1-i) * B;
end
A_x = [A_x; H_target * Phi_N];
A_U = [A_U; H_target * Gamma_N];
b_tot = [b_tot; h_target];

% Ora applichiamo la slice: x_0(slice_dims) = slice_vals
% x_0 = M_plot * x_plot + M_slice * slice_vals
M_plot = zeros(n, n_plot);
for idx = 1:n_plot
    M_plot(plot_dims(idx), idx) = 1;
end

M_slice = zeros(n, length(slice_dims));
for idx = 1:length(slice_dims)
    M_slice(slice_dims(idx), idx) = 1;
end

% Sostituiamo in A_x * x_0
% A_x * (M_plot * x_plot + M_slice * slice_vals) + A_U * U <= b_tot
% (A_x * M_plot) * x_plot + A_U * U <= b_tot - A_x * M_slice * slice_vals

A_plot = A_x * M_plot;
A_U_tot = A_U;
b_new = b_tot - A_x * M_slice * slice_vals;

% Il nuovo poliedro ha variabili Z = [x_plot; U]
A_poly = [A_plot, A_U_tot];
poly_full = Polyhedron(A_poly, b_new);

% Proiettiamo SOLTANTO sulle prime n_plot variabili (che sono x_plot, tipicamente 3D)
% Questa proiezione MPT3 da (3+N*m) a 3D è rapidissima rispetto a 6D
poly_full.minHRep();
poly_slice_out = projection(poly_full, 1:n_plot);
poly_slice_out.minHRep();

end
