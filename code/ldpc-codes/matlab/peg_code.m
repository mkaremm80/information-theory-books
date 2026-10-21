function code = peg_code(n, m, dv)
%PEG_CODE Progressive edge growth (Hu, Eleftheriou, Arnold 2005). dv: scalar or
% vector of variable degrees. Ties are broken with MATLAB's random stream.
dv = dv(:).' .* ones(1, n);
chk = cell(m, 1); var = cell(n, 1); cdeg = zeros(m, 1);
[~, order] = sort(dv);
for j = order
    for k = 1:dv(j)
        if k == 1
            cand = find(cdeg == min(cdeg));
        else
            seen = false(m, 1); seen(var{j}) = true;
            frontier = var{j};
            while true
                nv = unique([chk{frontier}]);
                nc = setdiff(unique([var{nv}]), find(seen));
                if isempty(nc) || all(seen)
                    cand = find(~seen); break;
                end
                seen(nc) = true; frontier = nc;
                if all(seen)
                    cand = nc(:); break;
                end
            end
            cand = cand(cdeg(cand) == min(cdeg(cand)));
        end
        c = cand(randi(numel(cand)));
        chk{c}(end + 1) = j; var{j}(end + 1) = c; cdeg(c) = cdeg(c) + 1;
    end
end
cn = [var{:}].'; vn = repelem((1:n).', cellfun(@numel, var));
code = ldpc_from_edges(m, n, cn, vn, sprintf('PEG(%d,%d)', n, n - m));
end
