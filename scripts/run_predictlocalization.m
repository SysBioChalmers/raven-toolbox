% run_predictlocalization.m -- run RAVEN predictLocalization N times on yeast-GEM + DeepLoc GSS.
%
% Head-to-head baseline for the deterministic MILP. predictLocalization is stochastic (random moves,
% no rng seed, wall-clock budget, greedy acceptance), so we run it N times to characterise its
% run-to-run variance. Reads the GSS table written by `compare_predictlocalization.py --prep` and
% writes one gene->compartment CSV per run plus a timing/score meta CSV.
%
% Override before calling (matlab -batch "N=5; maxT=2; run('scripts/run_predictlocalization.m')"):
%   ravenDir, modelFile, gssFile, outDir, N, maxT, removeExch
if ~exist('ravenDir','var');  ravenDir  = 'C:/Work/GitHub/RAVEN'; end
if ~exist('modelFile','var'); modelFile = 'C:/Work/GitHub/yeast-GEM/model/yeast-GEM.yml'; end
if ~exist('gssFile','var');   gssFile   = '.research_tmp/pl/deeploc_gss_yeast.csv'; end
if ~exist('outDir','var');    outDir    = '.research_tmp/pl'; end
if ~exist('N','var');         N    = 5; end
if ~exist('maxT','var');      maxT = 2; end
if ~exist('removeExch','var'); removeExch = true; end
% Transport cost per metabolite. predictLocalization's default (0.5) is sized for clean
% integer-style scores and overwhelms soft DeepLoc probabilities (every gene stays in `c`); dial it
% down into the per-metabolite-per-compartment range of the score deltas, as the yeast localization
% benchmark recommends. Matched to the MILP arm's transport_cost for a fair head-to-head.
if ~exist('transCost','var'); transCost = 0.05; end

addpath(genpath(ravenDir));

% --- load model (YAML = pure MATLAB, no libSBML); fall back to SBML import ---
try
    model = readYAMLmodel(modelFile);
catch
    model = importModel(strrep(modelFile, '.yml', '.xml'));
end

% predictLocalization documents that exchange/demand/sink reactions should be removed first.
if removeExch
    [~, exchIdx] = getExchangeRxns(model);
    model = removeReactions(model, exchIdx, true, true, true);
end
fprintf('model: %d rxns, %d genes, %d comps\n', numel(model.rxns), numel(model.genes), numel(model.comps));

% --- build GSS from the CSV (gene_id + one column per compartment; scores already normalised) ---
T = readtable(gssFile, 'ReadVariableNames', true, 'Delimiter', ',');
GSS.genes        = cellstr(string(T{:,1}));
GSS.compartments = T.Properties.VariableNames(2:end)';
GSS.scores       = T{:,2:end};
fprintf('GSS: %d genes x %d compartments (%s)\n', numel(GSS.genes), numel(GSS.compartments), strjoin(GSS.compartments,','));

meta = [];
for i = 1:N
    % predictLocalization can crash in its epilogue if the wall-clock loop ends mid-iteration
    % (an undefined geneScore at line ~590) -- a real robustness bug. Guard each run so one bad
    % draw does not kill the whole batch; geneLoc is written per-run so partials survive.
    t0 = tic;
    try
        [~, geneLoc, ~, sc] = predictLocalization(model, GSS, 'c', 'maxTime', maxT, ...
                                                  'transportCost', transCost);
    catch ME
        fprintf('run %d FAILED (%s) after %.0fs -- skipped\n', i, ME.message, toc(t0));
        continue;
    end
    secs = toc(t0);
    Tg = table(geneLoc.genes(:), geneLoc.comps(:), 'VariableNames', {'gene','comp'});
    writetable(Tg, fullfile(outDir, sprintf('geneloc_run_%03d.csv', i)));
    meta = [meta; i, maxT, secs, sc.totScore, sc.geneScore, sc.transCost]; %#ok<AGROW>
    fprintf('run %d: %.0fs, %d genes localized, totScore=%.1f geneScore=%.1f transCost=%.1f\n', ...
            i, secs, numel(geneLoc.genes), sc.totScore, sc.geneScore, sc.transCost);
end
Tm = array2table(meta, 'VariableNames', {'run','maxTime_min','seconds','totScore','geneScore','transCost'});
writetable(Tm, fullfile(outDir, 'pl_runs_meta.csv'));
fprintf('ALL DONE (%d runs)\n', N);
