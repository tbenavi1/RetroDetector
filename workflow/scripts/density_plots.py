from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

tp_geneids = set()
with open(snakemake.input[0], "r") as input_tp_file:
	for line in input_tp_file:
		geneid = line.strip()
		tp_geneids.add(geneid)

fn_geneids = set()
with open(snakemake.input[1], "r") as input_fn_file:
	for line in input_fn_file:
		geneid = line.strip()
		fn_geneids.add(geneid)

geneid_to_idents = defaultdict(list)
with open(snakemake.input[2], "r") as input_retrogene_file:
	for line in input_retrogene_file:
		geneid = line.strip().split()[0]
		pct_ident = float(line.strip().split()[7])
		geneid_to_idents[geneid].append(pct_ident)

tp_idents = []
fn_idents = []
for geneid in geneid_to_idents:
	idents = geneid_to_idents[geneid]
	avg_ident = sum(idents)/len(idents)
	if geneid in tp_geneids:
		tp_idents.append(avg_ident)
	else:
		assert geneid in fn_geneids
		fn_idents.append(avg_ident)

tp_density = gaussian_kde(tp_idents)
fn_density = gaussian_kde(fn_idents)
x_vals = range(0, 101, 1)
plt.plot(x_vals, tp_density(x_vals), label="Long and short reads")
plt.plot(x_vals, fn_density(x_vals), label="Long reads only")
plt.xlabel("Retrogene percent identity")
plt.ylabel("Density")
plt.legend()
plt.savefig(snakemake.output[0])
plt.close()
