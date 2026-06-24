function Q = gram_schmidt(V)
    [n, m] = size(V);
    Q = zeros(n, m);
    k = 0;
    for j = 1:m
        v = V(:, j);
        for i = 1:k
            v = v - (Q(:,i)' * v) * Q(:,i);
        end
        if norm(v) > 1e-10
            k = k + 1;
            Q(:, k) = v / norm(v);
        end
    end
    Q = Q(:, 1:k);
end