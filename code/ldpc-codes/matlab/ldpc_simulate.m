function res = ldpc_simulate(encode, decode, k, R, ebn0_list, info, max_frames, min_errors, batch)
%LDPC_SIMULATE BER/FER versus Eb/N0. encode(u) -> [codeword, sent bits];
% decode(llr of sent bits) -> [decoded codeword, iterations]; info = positions of
% the message bits in the decoded codeword (1-based).
res.ebn0 = ebn0_list; res.ber = []; res.fer = []; res.frames = []; res.avg_iter = [];
for eb = ebn0_list
    be = 0; fe = 0; fr = 0; its = 0;
    while fr < max_frames && fe < min_errors
        u = randi([0 1], batch, k);
        [~, x] = encode(u);
        [ch, it] = decode(awgn_llr(x, eb, R));
        e = sum(ch(:, info) ~= u, 2);
        be = be + sum(e); fe = fe + sum(e > 0); fr = fr + batch; its = its + sum(it);
    end
    res.ber(end + 1) = be / (fr * k); res.fer(end + 1) = fe / fr;
    res.frames(end + 1) = fr; res.avg_iter(end + 1) = its / fr;
    fprintf('Eb/N0 = %4.2f dB  FER = %.2e  BER = %.2e  frames = %d  avg it = %.1f\n', ...
            eb, res.fer(end), res.ber(end), fr, res.avg_iter(end));
end
end
