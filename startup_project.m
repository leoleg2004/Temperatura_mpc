function startup_project()
    % Aggiunge tutte le sottocartelle del progetto al path di MATLAB
    project_dir = fileparts(mfilename('fullpath'));
    addpath(genpath(project_dir));
    
    % FIX PER I TOOLBOX MPT3/YALMIP (risolve l'errore delle cartelle spostate)
    toolboxes_dir = '/Users/leonardoleggeri/Desktop/PROGETTI/Automazione/Multivaribaile/Esercitazioni/toolboxes';
    if isfolder(toolboxes_dir)
        addpath(genpath(toolboxes_dir));
        % Avvia MPT3 silenziosamente se non è già attivo
        if exist('mpt_init', 'file') == 2
            evalc('mpt_init'); 
        end
    else
        disp('ATTENZIONE: Cartella toolboxes MPT3 non trovata nel percorso specificato!');
    end
    
    disp('Path inizializzato correttamente per Progetto_MPC_Calore_Tre_Camere (incluso MPT3).');
end
