from collections import defaultdict

long_truepositives = set()
with open(snakemake.input[0], "r") as input_truepositives_file:
	for line in input_truepositives_file:
		geneid = line.strip()
		long_truepositives.add(geneid)

all_spanned = set()
singly_spanned = set()
multiply_spanned = set()
geneid_to_num_spanned = defaultdict(int)
with open(snakemake.input[1], "r") as input_spannedjunctions_file, open(snakemake.output[0], "w") as output_file:
	for line in input_spannedjunctions_file:
		geneid, transcript, num_junctions, num_spanned_junctions = line.strip().split()
		num_junctions, num_spanned_junctions = int(num_junctions), int(num_spanned_junctions)
		previous_num_spanned_junctions = geneid_to_num_spanned[geneid]
		geneid_to_num_spanned[geneid] = max(previous_num_spanned_junctions, num_spanned_junctions)
	for geneid in geneid_to_num_spanned:
		all_spanned.add(geneid)
		num_spanned = geneid_to_num_spanned[geneid]
		if num_spanned > 1:
			multiply_spanned.add(geneid)
		else:
			singly_spanned.add(geneid)
	output_file.write(f"support\tspanned\t# True Positives\t#False Positives\t#False Negatives\n")
	true_positives = len(all_spanned.intersection(long_truepositives))
	false_positives = len(all_spanned - long_truepositives)
	false_negatives = len(long_truepositives - all_spanned)
	#print(long_truepositives - all_spanned)
	#print(false_negatives)
	output_file.write(f"all support\tsingly or multiply spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")
	true_positives = len(singly_spanned.intersection(long_truepositives))
	false_positives = len(singly_spanned - long_truepositives)
	false_negatives = len(long_truepositives - singly_spanned)
	output_file.write(f"all support\tsingly spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")
	true_positives = len(multiply_spanned.intersection(long_truepositives))
	false_positives = len(multiply_spanned - long_truepositives)
	false_negatives = len(long_truepositives - multiply_spanned)
	output_file.write(f"all support\tmultiply spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")


no_non_overlapping_all_spanned = set()
no_non_overlapping_singly_spanned = set()
no_non_overlapping_multiply_spanned = set()
no_non_overlapping_geneid_to_num_spanned = defaultdict(int)
with open(snakemake.input[2], "r") as input_no_non_overlapping_spannedjunctions_file, open(snakemake.output[0], "a") as output_file, open(snakemake.output[1], "w") as output_tp_file, open(snakemake.output[2], "w") as output_fn_file:
	for line in input_no_non_overlapping_spannedjunctions_file:
		geneid, transcript, num_junctions, num_spanned_junctions = line.strip().split()
		num_junctions, num_spanned_junctions = int(num_junctions), int(num_spanned_junctions)
		previous_num_spanned_junctions = no_non_overlapping_geneid_to_num_spanned[geneid]
		no_non_overlapping_geneid_to_num_spanned[geneid] = max(previous_num_spanned_junctions, num_spanned_junctions)
	for geneid in no_non_overlapping_geneid_to_num_spanned:
		no_non_overlapping_all_spanned.add(geneid)
		num_spanned = no_non_overlapping_geneid_to_num_spanned[geneid]
		if num_spanned > 1:
			no_non_overlapping_multiply_spanned.add(geneid)
		else:
			no_non_overlapping_singly_spanned.add(geneid)
	true_positives = len(no_non_overlapping_all_spanned.intersection(long_truepositives))
	for geneid in no_non_overlapping_all_spanned.intersection(long_truepositives):
		output_tp_file.write(f"{geneid}\n")
	false_positives = len(no_non_overlapping_all_spanned - long_truepositives)
	false_negatives = len(long_truepositives - no_non_overlapping_all_spanned)
	for geneid in long_truepositives - no_non_overlapping_all_spanned:
		output_fn_file.write(f"{geneid}\n")
	output_file.write(f"no non-overlapping support\tsingly or multiply spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")
	true_positives = len(no_non_overlapping_singly_spanned.intersection(long_truepositives))
	#print(no_non_overlapping_singly_spanned.intersection(long_truepositives))
	false_positives = len(no_non_overlapping_singly_spanned - long_truepositives)
	false_negatives = len(long_truepositives - no_non_overlapping_singly_spanned)
	output_file.write(f"no non-overlapping support\tsingly spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")
	true_positives = len(no_non_overlapping_multiply_spanned.intersection(long_truepositives))
	print(no_non_overlapping_multiply_spanned.intersection(long_truepositives))
	false_positives = len(no_non_overlapping_multiply_spanned - long_truepositives)
	#print(no_non_overlapping_multiply_spanned - long_truepositives)
	false_negatives = len(long_truepositives - no_non_overlapping_multiply_spanned)
	output_file.write(f"no non-overlapping support\tmultiply spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")

no_alternate_all_spanned = set()
no_alternate_singly_spanned = set()
no_alternate_multiply_spanned = set()
no_alternate_geneid_to_num_spanned = defaultdict(int)
with open(snakemake.input[3], "r") as input_no_alternate_spannedjunctions_file, open(snakemake.output[0], "a") as output_file:
	for line in input_no_alternate_spannedjunctions_file:
		geneid, transcript, num_junctions, num_spanned_junctions = line.strip().split()
		num_junctions, num_spanned_junctions = int(num_junctions), int(num_spanned_junctions)
		previous_num_spanned_junctions = no_alternate_geneid_to_num_spanned[geneid]
		no_alternate_geneid_to_num_spanned[geneid] = max(previous_num_spanned_junctions, num_spanned_junctions)
	for geneid in no_alternate_geneid_to_num_spanned:
		no_alternate_all_spanned.add(geneid)
		num_spanned = no_alternate_geneid_to_num_spanned[geneid]
		if num_spanned > 1:
			no_alternate_multiply_spanned.add(geneid)
		else:
			no_alternate_singly_spanned.add(geneid)
	true_positives = len(no_alternate_all_spanned.intersection(long_truepositives))
	false_positives = len(no_alternate_all_spanned - long_truepositives)
	false_negatives = len(long_truepositives - no_alternate_all_spanned)
	output_file.write(f"no alternate support\tsingly or multiply spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")
	true_positives = len(no_alternate_singly_spanned.intersection(long_truepositives))
	false_positives = len(no_alternate_singly_spanned - long_truepositives)
	false_negatives = len(long_truepositives - no_alternate_singly_spanned)
	output_file.write(f"no alternate support\tsingly spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")
	true_positives = len(no_alternate_multiply_spanned.intersection(long_truepositives))
	false_positives = len(no_alternate_multiply_spanned - long_truepositives)
	#print(no_alternate_multiply_spanned - long_truepositives)
	false_negatives = len(long_truepositives - no_alternate_multiply_spanned)
	output_file.write(f"no alternate support\tmultiply spanned\t{true_positives}\t{false_positives}\t{false_negatives}\n")
