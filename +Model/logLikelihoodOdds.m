function llo = logLikelihoodOdds(params, data)

netSig = sqrt(params.var_x + params.var_s);

probPos = mog.create([+1, -1], [netSig, netSig], [params.p_match, 1-params.p_match]);
probNeg = mog.create([+1, -1], [netSig, netSig], [1-params.p_match, params.p_match]);

llo = mog.logpdf(data(:), probPos) - mog.logpdf(data(:), probNeg);
llo = reshape(llo, size(data));

end