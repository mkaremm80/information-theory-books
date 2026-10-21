%LAB11_ANOMALY  KL- and log-loss-based anomaly detection (MATLAB version).
%   Actions 1 = login, 2 = read, 3 = write, 4 = logout.
Pn = [0.05 0.80 0.10 0.05; 0.02 0.60 0.25 0.13; 0.02 0.50 0.30 0.18; 0.90 0.04 0.03 0.03];
Pa = [0.60 0.10 0.05 0.25; 0.05 0.30 0.55 0.10; 0.05 0.20 0.65 0.10; 0.85 0.05 0.05 0.05];
W = 50; rng(7);
train = markov_simulate(Pn, 20000, 1);
Ph = markov_estimate(train, 4, 1);              % smoothed ML baseline
q = markov_stationary(Ph);
freqw = @(w) (accumarray(w(:), 1, [4 1])' + 0.5) / (numel(w) + 2);
klsc  = @(w) kl_bits(freqw(w), q);
llsc  = @(w) -mean(log2(Ph(sub2ind([4 4], w(1:end-1), w(2:end)))));
win = @(s) reshape(s(1:floor(numel(s)/W)*W), W, [])';
valid = win(markov_simulate(Pn, 200*W, 1));
tn = win(markov_simulate(Pn, 200*W, 1));
ta = win(markov_simulate(Pa, 200*W, 1));
scores = {klsc, llsc}; names = {'KL', 'log-loss'};
fprintf('baseline action frequencies: %s\n', mat2str(round(q*1e3)/1e3));
for s = 1:2
    f = scores{s};
    v = sort(arrayfun(@(i) f(valid(i, :)), 1:size(valid, 1)));
    thr = v(min(numel(v), floor(0.99*numel(v)) + 1));
    fp = mean(arrayfun(@(i) f(tn(i, :)), 1:size(tn, 1)) > thr);
    fn = mean(arrayfun(@(i) f(ta(i, :)), 1:size(ta, 1)) <= thr);
    fprintf('%-8s threshold %.3f  false positive rate %.3f  false negative rate %.3f\n', names{s}, thr, fp, fn);
end
