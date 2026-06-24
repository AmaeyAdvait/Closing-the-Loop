function [T, n_v_co, n_v_c_obar, n_v_cbar_o, n_v_cbar_obar] = Kalman_Decomposition(Cm, Om)

    function C = my_setdiff(A, B)

    C = [];

        for i = 1:length(A)
            found = 0;
    
            for j = 1:length(B)
                if A(i) == B(j)
                    found = 1;
                    break;
                end
            end
    
            if ~found
                C = [C, A(i)];
            end
        end
    end

    function [R, pivot_cols] = rref_manual(A, tol)
        if nargin < 2
            tol = 1e-6; 
        end

        R = double(A);
        [m, n] = size(R);
        pivot_cols = [];
        row = 1;
        for col = 1:n
            if row > m
                break; 
            end
            [mx, idx] = max(abs(R(row:m, col)));
            if mx < tol
                continue; 
            end
            idx = idx + row - 1;
            pivot_cols(end+1) = col;
            R([row, idx], :) = R([idx, row], :);
            R(row, :) = R(row, :) / R(row, col);
            for r = 1:m
                if r ~= row
                    R(r, :) = R(r, :) - R(r, col) * R(row, :);
                end
            end
            row = row + 1;
        end
    end

    function N = null_manual(A, tol)
        if nargin < 2
            tol = 1e-6; 
        end

        [~, n] = size(A);
        [R, pc] = rref_manual(A, tol);
        free_cols = my_setdiff(1:n, pc);
        N = zeros(n, length(free_cols));
        for k = 1:length(free_cols)
            fc = free_cols(k);
            vec = zeros(n, 1);
            vec(fc) = 1;
            for j = 1:length(pc)
                vec(pc(j)) = -R(j, fc);
            end
            N(:, k) = vec;
        end
    end

    function B = col_space(A, tol)
        if nargin < 2
            tol = 1e-6; 
        end

        [~, pc] = rref_manual(A, tol);
        if isempty(pc)
            B = zeros(size(A,1), 0);
        else
            B = A(:, pc);
        end
    end

    function Z_I = Zassenhaus_Intersection(U, W)

        [u1,u2] = size(U);
        [~,w2] = size(W);

        if u2 == 0 || w2 == 0
            Z_I = zeros(u1, 0); 
            return;
        end

        S = [U',   U';
             W',   zeros(w2, u1)]; %Constructed Matrix
        R = rref_manual(S, 1e-6);

        left_half  = R(:, 1:u1);
        right_half = R(:, u1+1:2*u1);
        Z_I = zeros(u1, 0);
        for i = 1:size(R, 1)
            if norm(left_half(i,:)) < 1e-6 && norm(right_half(i,:)) > 1e-6
                Z_I = [Z_I, right_half(i,:)'];
            end
        end
    end

    % Four fundamental subspaces
    Im_C = col_space(Cm);       % controllable
    Null_CT = null_manual(Cm');    % uncontrollable
    Im_OT = col_space(Om');      % observable
    Null_O = null_manual(Om);     % unobservable

    fprintf('dim Im_C = %d\n', size(Im_C, 2));
    fprintf('dim Null_CT = %d\n', size(Null_CT, 2));
    fprintf('dim Im_OT = %d\n', size(Im_OT, 2));
    fprintf('dim Null_O = %d\n', size(Null_O, 2));

    % Four Kalman subspaces
    s_co = Zassenhaus_Intersection(Im_C, Im_OT);
    s_c_obar = Zassenhaus_Intersection(Im_C, Null_O);
    s_cbar_o = Zassenhaus_Intersection(Null_CT, Im_OT);
    s_cbar_obar = Zassenhaus_Intersection(Null_CT, Null_O);

    fprintf('dim s_co = %d\n', size(s_co, 2));
    fprintf('dim s_c_obar = %d\n', size(s_c_obar, 2));
    fprintf('dim s_cbar_o = %d\n', size(s_cbar_o, 2));
    fprintf('dim s_cbar_obar = %d\n', size(s_cbar_obar, 2));
    fprintf('sum = %d (should be %d)\n', ...
        size(s_co,2)+size(s_c_obar,2)+size(s_cbar_o,2)+size(s_cbar_obar,2), ...
        size(Cm,1));

    % Orthonormal bases
    v_co = gram_schmidt(s_co);
    v_c_obar = gram_schmidt(s_c_obar);
    v_cbar_o = gram_schmidt(s_cbar_o);
    v_cbar_obar = gram_schmidt(s_cbar_obar);

    n_v_co = size(v_co,2);
    n_v_c_obar = size(v_c_obar, 2);
    n_v_cbar_o = size(v_cbar_o, 2);
    n_v_cbar_obar = size(v_cbar_obar, 2);

    % Similarity Transformation T
    T = [v_c_obar, v_co, v_cbar_obar, v_cbar_o];

    [u1, cols] = size(T);
    fprintf('T: %d x %d\n', u1, cols);
end