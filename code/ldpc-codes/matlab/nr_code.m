function nr = nr_code(K, E, bg)
%NR_CODE 5G NR LDPC code block: lifted base graph (only the rows needed for the
% E transmitted bits), encoder and rate matching (rv = 0, no interleaving).
if nargin < 3, bg = []; end
[nr.bg, Z, ils] = nr_select(K, K / E, bg);
nr.K = K; nr.E = E; nr.Z = Z; nr.ils = ils;
V = nr_base_graph(nr.bg, ils);
P = V; P(V >= 0) = mod(V(V >= 0), Z);                  % shift values P = V mod Z
if nr.bg == 1, nr.kb = 22; else, nr.kb = 10; end
nr.Kp = nr.kb * Z;                                     % encoder input incl. fillers
nr.filler = K + 1:nr.Kp;
nfull = size(V, 2) * Z;
tx = setdiff(2 * Z + 1:nfull, nr.filler);              % first 2Z bits punctured
nr.Ncb = numel(tx);
if E <= numel(tx), tx = tx(1:E); end
ncol = max(nr.kb + 4, ceil(max(tx) / Z));              % base columns in use
nr.code = qc_lift(P(1:ncol - nr.kb, 1:ncol), Z);
nr.code.name = sprintf('NR BG%d K=%d E=%d Z=%d', nr.bg, K, E, Z);
nr.n = ncol * Z;
nr.enc = dual_encoder(nr.code, nr.Kp, 4 * Z);
nr.tx_pos0 = tx;
nr.tx_pos = tx(mod(0:E - 1, numel(tx)) + 1);           % circular buffer
end
