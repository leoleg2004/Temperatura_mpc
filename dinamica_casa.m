function dx = dinamica_casa(x, u)
% DINAMICA_CASA - Modello continuo non lineare della temperatura
%                 di un appartamento a tre stanze.
%
% Da utilizzare in un blocco "MATLAB Function" in Simulink.
%
% INPUTS:
%   x : vettore di stato [T1; T2; T3; Q1; Q2; Q3]
%   u : vettore degli ingressi (potenze di riferimento) [Q1_r; Q2_r; Q3_r]
%
% OUTPUTS:
%   dx: derivata dello stato rispetto al tempo [dT1; dT2; dT3; dQ1; dQ2; dQ3]

    % 1. Estrazione degli stati (Temperature e Potenze erogate)
    T1 = x(1);
    T2 = x(2);
    T3 = x(3);
    Q1 = x(4);
    Q2 = x(5);
    Q3 = x(6);

    % 2. Estrazione degli ingressi (Potenze di riferimento)
    Q1_r = u(1);
    Q2_r = u(2);
    Q3_r = u(3);

    % 3. Parametri del sistema
    % Capacità termiche [J/s]
    C1 = 6300;
    C2 = 4600;
    C3 = 4200;

    % Parametri di scambio termico di base [W/K]
    k1_bar = 16;
    k2_bar = 18;
    k3_bar = 19;
    kext = 9;   % Scambio con l'esterno

    % Costanti di tempo dei termosifoni [s]
    tau1 = 580;
    tau2 = 520;
    tau3 = 540;

    % Temperatura esterna costante [K]
    Text = 278;

    % 4. Calcolo dei coefficienti non lineari di scambio termico k_ij(t)
    % La formula dal progetto è: k_ij(t) = k_ij_bar + 4 / (1 + exp(-0.5 * ||Ti(t) - Tj(t)||_2))
    % Essendo scalari, la norma 2 è equivalente al valore assoluto.
    k12 = k1_bar + 4 / (1 + exp(-0.5 * abs(T1 - T2)));
    k13 = k2_bar + 4 / (1 + exp(-0.5 * abs(T1 - T3)));
    k23 = k3_bar + 4 / (1 + exp(-0.5 * abs(T2 - T3)));

    % 5. Equazioni differenziali per le temperature (dinamica delle stanze)
    dT1 = (Q1 - k12*(T1 - T2) - k13*(T1 - T3) - kext*(T1 - Text)) / C1;
    dT2 = (Q2 + k12*(T1 - T2) - k23*(T2 - T3) - kext*(T2 - Text)) / C2;
    dT3 = (Q3 + k13*(T1 - T3) + k23*(T2 - T3) - kext*(T3 - Text)) / C3;

    % 6. Equazioni differenziali per i termosifoni (dinamica degli attuatori)
    dQ1 = (Q1_r - Q1) / tau1;
    dQ2 = (Q2_r - Q2) / tau2;
    dQ3 = (Q3_r - Q3) / tau3;

    % 7. Composizione del vettore derivata
    dx = [dT1; dT2; dT3; dQ1; dQ2; dQ3];
end
