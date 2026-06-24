%Reduced Order Luenberger Observer

%Programmer's Log: This was the hardest block of code for me to personally
% implement, since the approach to find the L matrix was not explicitly
% discussed in class. Most of my time was spent trying to debug why L was
% a matrix of dimension 0x5, or why A12 was making Q negative semi-definite
% instead of positive. Very minor errors, but a lot of time went into
% debugging them.

function L = Luenberger_Observer(A12, A22)
    
    n = size(A22, 1);
    Q = A12' * A12;  % Q comes from the output coupling matrix

    % A22 must be stable for Gramian to converge
    ev = eig(A22);
    disp(ev);
    alpha = max(real(ev)) + 1;
    F = A22 - alpha.*eye(n);   % shift to make stable

    % Observability Gramian: solving F'P + PF = -Q
    % W_o = integral from 0 to inf of e^{F't} * Q * e^{Ft} dt
    ev_F = eig(F);
    T_end = -6 / max(real(ev_F));
    N = 1000;
    if mod(N,2) ~= 0
        N = N+1; 
    end
    dt = T_end / N;

    expF_dt  = expm(F.*dt);
    expFT_dt = expm(F'.*dt);
    P = zeros(n);
    expF_t = eye(n);
    expFT_t = eye(n);

    %Simpson's integral formula to compute the Observability Gramian

    for k = 0:N
        Wo = expFT_t * Q * expF_t;   % e^{F't}*Q*e^{Ft}: observability form

        if k == 0 || k == N
            w = 1;
        elseif mod(k,2) == 1
            w = 4;
        else
            w = 2;
        end

        P = P + w*Wo;

        expF_t  = expF_dt*expF_t;
        expFT_t = expFT_dt*expFT_t;
    end

    P = (dt/3) * P;
    P = (P + P') / 2;

    L = (P\A12');

    % Verify observer stability
    ev_obs = eig(A22 - L*A12); 
    % disp(ev_obs);
    if all(real(ev_obs) < 0)
        fprintf('Observer STABLE.\n');
    else
        fprintf('Unstable - increase alpha.\n');
    end

    L = L';
end