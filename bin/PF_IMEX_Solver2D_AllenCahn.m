function [phi] = PF_IMEX_Solver2D_AllenCahn(LK,dx,dy,nx,ny,dt_phy,phi,s,mask)

%Enforce mask to be logical
mask        = logical(mask);

%Number of phases
Np          = size(phi,3);

%LK vectorize
LKvec_full  = LK(:);

%Choose direct solver if the DOF is less than tot_direct
tol_direct  = 40000; 

%Loop through phases
for p = 1:Np
    %Current mask
    mask_p  = mask(:,:,p);
    %Find index of the mask_p > 0
    act_lin = find(mask_p);
    %Number of variable to be solved
    Nact    = numel(act_lin);
    %If there is no such phase, then go ahead
    if Nact == 0
        phi(:,:,p) = 0;continue;
    end

    % Get the row and column of act_lin
    [ia,ja] = ind2sub([ny,nx], act_lin); 

    % Choose the interior index for L/R/U/D
    hasL    = (ja > 1 );
    hasR    = (ja < nx);
    hasU    = (ia > 1 );
    hasD    = (ia < ny);

    % Generate the index of L/R/U/D
    indL  = zeros(Nact,1);  indL(hasL) = sub2ind([ny,nx], ia(hasL)  , ja(hasL)-1);
    indR  = zeros(Nact,1);  indR(hasR) = sub2ind([ny,nx], ia(hasR)  , ja(hasR)+1);
    indU  = zeros(Nact,1);  indU(hasU) = sub2ind([ny,nx], ia(hasU)-1, ja(hasU)  );
    indD  = zeros(Nact,1);  indD(hasD) = sub2ind([ny,nx], ia(hasD)+1, ja(hasD)  );

    % Calculate coefficients for C/L/R/U/D
    LKc   = LKvec_full(act_lin);
    coefL = zeros(Nact,1); coefL(hasL) = -(LKvec_full(indL(hasL)) + LKc(hasL))/2 / dx^2;
    coefR = zeros(Nact,1); coefR(hasR) = -(LKvec_full(indR(hasR)) + LKc(hasR))/2 / dx^2;
    coefU = zeros(Nact,1); coefU(hasU) = -(LKvec_full(indU(hasU)) + LKc(hasU))/2 / dy^2;
    coefD = zeros(Nact,1); coefD(hasD) = -(LKvec_full(indD(hasD)) + LKc(hasD))/2 / dy^2;
    coefC =-(coefL + coefR + coefU + coefD) + 1/dt_phy;
    
    % Build grid to active relation matrix
    g2c_local          = zeros(ny*nx,1,'int32');
    g2c_local(act_lin) = int32(1:Nact);   % compact index 从 1..Nact

    % Judge which is interior (int) and which is boundary (bnd)
    interiorFlag = (ia>1 & ia<ny & ja>1 & ja<nx);
    idx_int      = find( interiorFlag); 
    idx_bnd      = find(~interiorFlag);

    % Estimate the size of the rows
    est_nnz      = Nact + 4*numel(idx_int) + 2*numel(idx_bnd);
    rows         = zeros(est_nnz,1,'int32'); cols = zeros(est_nnz,1,'int32'); vals = zeros(est_nnz,1);
    
    %Initialize counter
    ptr          = 1;

    % Right hand size vector
    RHS          = zeros(Nact,1);

    % ---------------------------
    % CREATE INTERIOR MATRIX'S L
    % ---------------------------
    if ~isempty(idx_int)
        %Get the ID of interior
        compIdxAll = (1:Nact)';
        compIdxInt = compIdxAll(idx_int);
        %Diagonal
        rows(ptr:ptr+numel(compIdxInt)-1) = compIdxInt;
        cols(ptr:ptr+numel(compIdxInt)-1) = compIdxInt;
        vals(ptr:ptr+numel(compIdxInt)-1) = coefC(idx_int);
        ptr                               = ptr + numel(compIdxInt);

        % Left neighbours
        neighL_comp = g2c_local(indL(idx_int));
        condL       = (hasL(idx_int) & (neighL_comp > 0));
        nL          = nnz(condL);
        if nL > 0
            rows(ptr:ptr+nL-1) = compIdxInt(condL);
            cols(ptr:ptr+nL-1) = double(neighL_comp(condL));
            vals(ptr:ptr+nL-1) = coefL(idx_int(condL));
            ptr = ptr + nL;
        end

        % Right neighbours
        neighR_comp = g2c_local(indR(idx_int));
        condR       = (hasR(idx_int) & (neighR_comp > 0));
        nR          = nnz(condR);
        if nR > 0
            rows(ptr:ptr+nR-1) = compIdxInt(condR);
            cols(ptr:ptr+nR-1) = double(neighR_comp(condR));
            vals(ptr:ptr+nR-1) = coefR(idx_int(condR));
            ptr = ptr + nR;
        end

        % Upper neighbours
        neighU_comp = g2c_local(indU(idx_int));
        condU       = (hasU(idx_int) & (neighU_comp > 0));
        nU          = nnz(condU);
        if nU > 0
            rows(ptr:ptr+nU-1) = compIdxInt(condU);
            cols(ptr:ptr+nU-1) = double(neighU_comp(condU));
            vals(ptr:ptr+nU-1) = coefU(idx_int(condU));
            ptr = ptr + nU;
        end

        % Lower neighbours
        neighD_comp = g2c_local(indD(idx_int));
        condD       = (hasD(idx_int) & (neighD_comp > 0));
        nD          = nnz(condD);
        if nD > 0
            rows(ptr:ptr+nD-1) = compIdxInt(condD);
            cols(ptr:ptr+nD-1) = double(neighD_comp(condD));
            vals(ptr:ptr+nD-1) = coefD(idx_int(condD));
            ptr = ptr + nD;
        end

        % RHS
        phi_p_vec       = phi(:,:,p); phi_p_vec = phi_p_vec(:);
        s_p_vec         = s{p}(:);
        RHS(compIdxInt) = phi_p_vec(act_lin(idx_int))/dt_phy + s_p_vec(act_lin(idx_int));
    end

    % ---------------------------
    % ADD BOUNDARY
    % ---------------------------
    if ~isempty(idx_bnd)
        compIdxAll = (1:Nact)';
        for tt = 1:numel(idx_bnd)
            pos      = idx_bnd(tt); 
            comp_row = compIdxAll(pos);

            i = ia(pos); j = ja(pos);
            % Choose neighbours
            if i == 1
                ni = i+1; nj = j;
            elseif i == ny
                ni = i-1; nj = j;
            elseif j == 1
                ni = i; nj = j+1;
            elseif j == nx
                ni = i; nj = j-1;
            else
                ni = i; nj = j;
            end

            neigh_lin  = sub2ind([ny,nx], ni, nj);
            neigh_comp = g2c_local(neigh_lin);  

            % row:  phi(center) - phi(neigh_inside) = 0
            rows(ptr) = comp_row; cols(ptr) = comp_row; vals(ptr) = 1.0; ptr = ptr + 1;
            if neigh_comp > 0
                rows(ptr) = comp_row; cols(ptr) = double(neigh_comp); vals(ptr) = -1.0; ptr = ptr + 1;
            else
                % neighbor
            end
        end
    end

    % Trim triplet arrays
    rows = rows(1:ptr-1); cols = cols(1:ptr-1); vals = vals(1:ptr-1);

    % Build sparse matrix
    L    = sparse(double(rows), double(cols), vals, Nact, Nact);

    if Nact <= tol_direct
        x = L \ RHS;
    else
        % bicgstab solver
        try
            setup_opts.type = 'ilutp'; setup_opts.droptol = 1e-3;
            [Lpre,Upre]     =  ilu(L, setup_opts);
            [x,flag]        =  bicgstab(L, RHS, 1e-10, 500, Lpre, Upre);
            if flag ~= 0
                warning('bicgstab not converged，use direct solver');
                x = L \ RHS;
            end
        catch
            x = L \ RHS;
        end
    end

    % Reasign
    phi_p_vec          = zeros(ny*nx,1);
    phi_p_vec(act_lin) = x;
    tmp                = reshape(phi_p_vec, ny, nx);
    tmp(~mask_p)       = 0;
    phi(:,:,p)         = tmp;
end

end
