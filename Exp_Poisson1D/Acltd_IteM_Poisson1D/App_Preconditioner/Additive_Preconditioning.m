function z = Additive_Preconditioning(idx_sets, sub_precs, num_sub, coarse_prec, r)
    z = zeros(size(r,1),1);
    for j = 1:num_sub
        idx = idx_sets{j};
        rj = r(idx);
        zj = sub_precs{j}(rj);
        z(idx) = z(idx) + zj;
    end
    z = z + coarse_prec(r);
end