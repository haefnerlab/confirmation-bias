function s = lloToEvidence(params, llos)

gen_sig_s = sqrt(Model.getEvidenceVariance(params.sensory_info));
p_match = params.category_info;

num_grid = 200;
sgrid = sort([linspace(-1-5*gen_sig_s, -1+5*gen_sig_s, num_grid) linspace(+1-5*gen_sig_s, +1+5*gen_sig_s, num_grid)])';

mog_pos = mog.create([-1 +1], [gen_sig_s gen_sig_s], [1-p_match p_match]);
mog_neg = mog.create([-1 +1], [gen_sig_s gen_sig_s], [p_match 1-p_match]);

llogrid = mog.logpdf(sgrid, mog_pos) - mog.logpdf(sgrid, mog_neg);

s = interp1(llogrid, sgrid, llos);

s = reshape(s, size(llos));

end