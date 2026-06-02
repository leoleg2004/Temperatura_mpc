function [Ad, Bd] = discretizza_modello(Ac, Bc, nx, nu, Ts)
    % Converte il sistema da tempo continuo a discreto
    sys_c = ss(Ac, Bc, eye(nx), zeros(nx, nu));
    sys_d = c2d(sys_c, Ts, 'zoh'); % Zero-Order Hold
    Ad = sys_d.A;
    Bd = sys_d.B;
end
