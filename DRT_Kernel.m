function [K_adj,W_adj,D_adj] = DRT_Kernel(type_intercept,type_kernel,x,ppd_t,t,y_total)

    M =length(x); % length of the experimental data vectors
    dt=1/ppd_t;
    N=length(t);
    
    % K_matrix
    K = zeros(M,N);
    if type_kernel ==1  % Finite Warburg kernel (transmissive boundary)
        for m =1:M
            for n=1:N
                K(m,n) = dt*sqrt(1i*exp(t(n)-x(m)))*coth(sqrt(1i*exp(t(n)-x(m)))); 
            end
        end
    elseif type_kernel == 2 % Bounded Warburg kernel (blocking boundary)
        for m =1:M
            for n=1:N
                if (t(n)-x(m)) > 6
                % K_D(m,n) = dtD*sqrt(1i*exp(tD(n)-x(m)))*(tanh(sqrt(1i*exp(tD(n)-x(m))))); 
                K(m,n) = dt*((1/sqrt(2*exp(t(n)-x(m))))*(1-1i)-(1/(2*exp(t(n)-x(m))))*1i)^-1;
                
                else
                K(m,n) = dt*sqrt(1i*exp(t(n)-x(m)))*(besseli(1,(sqrt(1i*exp(t(n)-x(m)))))/besseli(0,(sqrt(1i*exp(t(n)-x(m))))));     
                end
            end
        end
    elseif type_kernel == 0 % RC (relaxation) kernel
        for m =1:M
            for n=1:N
                K(m,n) =dt*(1+1i*exp(t(n)-x(m)))^-1; % Z 
            end
        end
    end
    
    % W_matrix (weight) % It is one weighting for R and D.
    W = zeros(M);
    for m = 1:M
        %W(m,m) = 1; % only diagonal components
        W(m,m) = abs(y_total(m))^-1; % only diagonal components
    end
    
    % D_matrix
    
    D = zeros(N,N);
    for n = 1:N
        if n ==1 % forward
        D(n,n)=1;
        D(n,n+1)= -2;
        D(n,n+2)=1;
        elseif n == N % backward
        D(n,n-2)=1;
        D(n,n-1)= -2;
        D(n,n)=1;    
        else % central
        D(n,n-1)=1;
        D(n,n)= -2;
        D(n,n+1)=1;
        end
    end
    
    D = D/(dt);


    %% Adjust for adding an unknown resistance
    if type_intercept == 1
        K_adj = [K,ones(M,1)];
        W_adj = W;
        D_adj = zeros(N+1,N+1); D_adj(1:N,1:N)=D; 
    elseif type_intercept == 0
        K_adj = K;
        W_adj = W;
        D_adj = D;
    end

end