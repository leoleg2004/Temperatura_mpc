function [Fx, fx, Fu, fu, dU_min, dU_max, Gx, gx] = imposta_vincoli(X_min, X_max, U_min, U_max, Ts)
    % IMPOSTA_VINCOLI Costruisce le matrici poliedriche dei vincoli per l'MPC.
    % È una funzione generalizzata per qualsiasi sistema.
    %
    % Input:
    %   X_min, X_max: Limiti inferiori e superiori degli stati (nx x 1)
    %   U_min, U_max: Limiti inferiori e superiori degli ingressi (nu x 1)
    %   Rate_max: Massimo rateo di variazione degli ingressi (nu x 1)
    %   Ts: Tempo di campionamento
    %
    % Output:
    %   Fx, fx: Matrici vincoli di stato Fx * x <= fx
    %   Fu, fu: Matrici vincoli di ingresso Fu * u <= fu
    %   dU_min, dU_max: Vincoli di rateo
    %   Gx, gx: Copia di Fx, fx (per retrocompatibilità)

    nx = length(X_min);
    nu = length(U_min);

    % Vincoli di Rateo
    dU_max = U_max; 
    dU_min = U_min;

    % Vincoli Poliedrici di Stato
    Fx = [eye(nx); -eye(nx)]; 
    fx = [X_max; -X_min];

    % Vincoli Poliedrici di Ingresso
    Fu = [eye(nu); -eye(nu)]; 
    fu = [U_max; -U_min];
    
    % Per retrocompatibilità con altri codici MPC pre-esistenti
    Gx = Fx; 
    gx = fx; 
end
