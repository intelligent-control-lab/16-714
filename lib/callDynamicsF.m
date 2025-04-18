function out = callDynamicsF(f, x, u, t)
    if nargin(f) == 2
        out = f(x, u);
    else
        out = f(x, u, t);
    end
end