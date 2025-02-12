function B = bspline(x, xl, xr, ndx, bdeg)
    dx = (xr - xl) / ndx;
    t = xl + dx * [-bdeg:ndx-1];
    T = (0 * x + 1) * t;
    X = x * (0 * t + 1);
    P = (X - T) / dx;
    B = (T <= X) & (X < (T + dx));
    r = [2:length(t) 1];
    for k = 1:bdeg
        B = (P .* B + (k + 1 - P) .* B(:, r)) / k;
    end
end