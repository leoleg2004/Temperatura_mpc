function [G_inf, g_inf] = cis(A, B, x_ref, u_ref, Fx, fx, Fu, fu, Q, R)
    % Control Invariant Set Iterativo tramite Metodo Matriciale (No Polyhedron)
    [K, ~, ~] = dlqr(A, B, Q, R);
    A_cl = A - B*K;
    
    % Traslazione dei vincoli rispetto al riferimento z = x - x_ref
    Gx = Fx; 
    gx = fx - Fx*x_ref;
    
    Gu_x = [-K; K]; 
    gu_x = [fu(1:3) - Fu(1:3,:)*u_ref; fu(4:6) - Fu(4:6,:)*u_ref];
    
    G = [Gx; Gu_x]; 
    g = [gx; gu_x];

    G_inf = G; g_inf = g;
    max_iter = 100; tol = 1e-6;
    for i = 1:max_iter
        G_next = G * (A_cl^i);
        if norm(G_next, inf) < tol
            fprintf('CIS convergente! L''algoritmo ha trovato O_inf in %d iterazioni.\n', i);
            break;
        end
        G_inf = [G_inf; G_next];
        g_inf = [g_inf; g];
        if i == max_iter
            fprintf('Attenzione: il calcolo del CIS non è arrivato a convergenza esatta dopo %d iterazioni.\n', max_iter);
        end
    end
    
    % Ritraslazione di g_inf in modo che G_inf * x <= g_inf 
    % (dato che attualmente G_inf * z <= g_inf)
    g_inf = g_inf + G_inf * x_ref;
end
