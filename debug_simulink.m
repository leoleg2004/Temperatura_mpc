startup_project();
try
    mdl = 'Stanza';
    load_system(mdl);
    opspec = operspec(mdl);
    disp('Operspec created successfully.');
    
    % Try to compile the model to see the error
    [sys, x0, str, ts] = sys_c(mdl); 
catch e
    disp('ERROR DETECTED:');
    disp(e.message);
    if ~isempty(e.cause)
        for i=1:length(e.cause)
            disp(['Cause ', num2str(i), ': ', e.cause{i}.message]);
        end
    end
end
