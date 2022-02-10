%% Impulse Response Decomposition 
% Consumption

%% Iteration

%--------------------------------------------------------------------------
% Backward Induction (Value Function)
%--------------------------------------------------------------------------
Apath = zeros(Naz,Naz,N);
cpath = zeros(Na,Nz,N);
V = VSS;

C_case = zeros(N,2);

for i_case = 1:2
    
    for n = N:-1:1
        
        if i_case == 1
            r = TRAN.rpath(n);
            w = wSS;
        else
            w = TRAN.wpath(n);
            r = rSS;
        end

        %----------------------------------
        % Forward Difference 
        %----------------------------------
        dVaF(1:Na-1,:) = (V(2:Na,:) - V(1:Na-1,:))/da;
        dVaF(Na,:) = (r*amax + w*zgrid).^(-gamma);
        cF = dVaF.^(-1/gamma);
        sF = r*amatrix + w*zmatrix - cF;

        %----------------------------------
        % Backward Difference 
        %----------------------------------
        dVaB(2:Na,:) = (V(2:Na,:) - V(1:Na-1,:))/da;
        dVaB(1,:) = (r*amin + w*zgrid).^(-gamma);
        cB = dVaB.^(-1/gamma);
        sB = r*amatrix + w*zmatrix - cB;

        %----------------------------------
        % Upwind Scheme
        %----------------------------------
        Upwind_F = sF > 0;
        Upwind_B = sB < 0;
        Upwind0 = 1 - Upwind_F - Upwind_B;

        c0 = r*amatrix + w*zmatrix;
        dVa0 = c0.^(-gamma);

        c = Upwind_F.*cF + Upwind_B.*cB + Upwind0.*c0;
        cpath(:,:,n) = c;
        dV0 = c0.^(-gamma);
        s = r*amatrix + w*zmatrix - c;
        u = c.^(1-gamma)/(1-gamma);

        %----------------------------------
        % Construct A Matrix
        %----------------------------------
        A_up = max(sF,0)/da;
        A_cent = -max(sF,0)/da + min(sB,0)/da;
        A_low = -min(sB,0)/da;

        A_low(1,:) = 0;
        A_low = A_low(:);
        A_up(Na,:) = 0;

        A0 = spdiags(A_cent(:),0,Naz,Naz) + spdiags(A_low(2:end),-1,Naz,Naz) + spdiags([0;A_up(:)],1,Naz,Naz);
        A = A0 + Zswitch;
        Apath(:,:,n) = A;

        Vprevious = ((1/dt+rho)*speye(Naz,Naz) - A)\(u(:) + V(:)/dt);
        V = reshape(Vprevious,Na,Nz);

    end
    %% Iterate KFE Forward 

    Cpath = CSS*ones(N,1);
    gpath = zeros(Naz,N);
    gpath(:,1) = gSS(:);

    for n = 1:N-1
        Atran = Apath(:,:,n)';
        g = (speye(Naz,Naz)/dt - Atran)\(gpath(:,n)/dt);
        gpath(:,n+1) = g;
        Cpath(n) = reshape(cpath(:,:,n),1,Naz)*g(:)*da;
    end
    
    C_case(:,i_case) = Cpath;
end
