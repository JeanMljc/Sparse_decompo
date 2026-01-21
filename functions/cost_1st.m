function[C] = cost_1st(M,p)
    % Cost matrix for convolutive dictionary
    % for 1ST order penalty
    %----------------------------------------------------------
    % INPUT 
    % M : number of atoms in the dictionary
    % p : distance used for Wp cost 
    %----------------------------------------------------------
    % OUTPUT 
    % C : cost matrix

    x = (1:M); % x is the discrete support
    
    C = abs(repmat(x,M,1)-repmat(x',1,M));

    if p ~= 1
        C = C.^p;     
    end
end