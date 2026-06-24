%Linear Quadratic Regulator

function K = LQR(A_co, B_co, C_co, n_v_c_o)

    Q_tilde = eye(size(C_co,1));
    R_tilde = eye(size(B_co,2));
    
    Q = C_co'*Q_tilde*C_co;
    R = -B_co*(R_tilde\B_co');
    
    H = [A_co, R; -Q, -A_co'];
    
    [v, ev] = eig(H);
    ev = diag(ev);
    jw = ((real(ev)) > -1e-6);
    
    v_new = v(:,~jw);
    X1 = v_new(1:n_v_c_o, :);
    X2 = v_new(n_v_c_o+1:2*n_v_c_o, :);
    
    P = X2/X1;
    K = R_tilde\(B_co'*P);
end