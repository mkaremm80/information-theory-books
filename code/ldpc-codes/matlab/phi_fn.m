function y = phi_fn(x)
%PHI_FN Gallager's function phi(x) = -ln tanh(x/2); it is its own inverse.
x = min(max(x, 1e-12), 50);
y = -log(tanh(x / 2));
end
