function sub = Get_OT_mass_disp(P, thr)
    % Extract the index in marginales between which the mass moves
    %----------------------------------------------------------
    % INPUT 
    % P : set of N OT plan (M x M x N)
    % thr : threshold below which we consider there is not mass
    %----------------------------------------------------------
    % OUTPUT 
    % sub : index in marginales between which the mass moves
    
    % Get size and number of OT tensors
    M = size(P,1);
    N = size(P,3) + 1;

    % Assign zero value to small displacement 
    P_thr = P;
    P_thr(P_thr<thr)=0;
    
    % Get the index between which the mass moves
    sub = cell(1,N-1);
    
    for n=1:N-1
        [sub1,sub2] = ind2sub([M,M],find(P_thr(:,:,n)~=0)); 
        % extract non-zero elements in P_thr
        sub{n} = [sub1,sub2];
    end
end