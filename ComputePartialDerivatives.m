% Compute Partial Derivatives time path

T = 200;
N = 800;
timespan = 0 + (T-0)*linspace(0,1,N)'.^1;
dt = T/N;

VarC = 2;
dC_v = zeros(N,VarC);

dvar = [0.0005,0.01];
PathSS = [rSS,wSS];

%--------------------------------------------------------------------------
% Compute r partial 
%--------------------------------------------------------------------------
V = VSS;
for n_var = 1:VarC
    
    for n2 = 1:N
        
        n2
        for n = n2:-1:1
            
            if n == n2
                pathnow = PathSS;
                pathnow(n_var) = pathnow(n_var) + dvar(n_var);
            else
                pathnow = PathSS;
            end
            
            r = pathnow(1);
            w = pathnow(2);
            
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
            Cpath(:,:,n) = c;
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

            Vprevious = ((1/dt+rho)*speye(Naz,Naz) - A)\(u(:) + V(:)/dt);
            V = reshape(Vprevious,Na,Nz);
        end

        % Compute initial consumption
        gnew = (speye(Na*Nz,Na*Nz) - dt*A')\gSS(:);
        dC = c(:)'*gnew(:)*da - CSS;
        dC_v(n2,n_var) = dC/dvar(n_var);
    end
end