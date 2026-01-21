function Arrow = Get_arrow(sub_OT)
    
    % Build arrow for plot
    %----------------------------------------------------------
    % INPUT
    % sub_OT : position of each arrow 
    %----------------------------------------------------------
    % OUTPUT 
    % Arrow : on the 2D grid (M x N)
    
    Arrow = [];
    
    for n=1:length(sub_OT)
        su = sub_OT{n};
        nP = length(su);
        AR = zeros(nP,4);
        AR(:,1) = n;
        AR(:,2) = su(:,1);
        AR(:,3) = 1;
        AR(:,4) = su(:,2)-su(:,1);
        Arrow = [Arrow; AR];
    end
end