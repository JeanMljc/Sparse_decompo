function [tra_dico] = perfect_decompo(mu, tra)
    % Recover the sequence of atoms in the dictionary that best fit to 
    % one given continuous (off-the-grid) trajectory
    %
    %----------------------------------------------------------
    % INPUT 
    % 
    % mu : center of atoms in the dictionary A
    % tra : center of IP in the trajectory in Y
    %----------------------------------------------------------
    % OUTPUT 
    % 
    % tra_dico : Trajectory in atoms space that best fit the real 
    % trajectory
    %----------------------------------------------------------
    [N_tra,N] = size(tra);
    tra_dico = zeros(N,N_tra);
    for i=1:N_tra
        for n=1:N
            dif = abs(mu-tra(i,n));
            [~,ind] = min(dif);
            tra_dico(n,i) = ind;
        end
    end
end
