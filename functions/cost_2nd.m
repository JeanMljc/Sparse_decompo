function[C] = cost_2nd(M,p)
    % Cost matrix for convolutive dictionary
    % for 2ND order penalty
    %----------------------------------------------------------
    % INPUT 
    % M : number of atoms in the dictionary
    % p : distance used for Wp cost 
    %----------------------------------------------------------
    % OUTPUT 
    % C : cost matrix
    
    [sub1,sub2,sub3] = ind2sub([M M M],(1:M^3));

    C = reshape(abs(2*sub2 - sub1 - sub3),M,M,M);
    
    if p ~= 1
        C = C.^p;
    end
end