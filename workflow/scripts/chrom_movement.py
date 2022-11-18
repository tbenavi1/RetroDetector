a_to_a = 0
a_to_x = 0
x_to_a = 0
x_to_x = 0
total = 0
with open(snakemake.input[0], "r") as input_file:
	for line in input_file:
		parent_location, retro_location = line.strip().split()[1:3]
		parent_chrom = parent_location.split(":")[0]
		retro_chrom = retro_location.split(":")[0]
		total += 1
		if parent_chrom != "X" and retro_chrom != "X":
			a_to_a += 1
		elif parent_chrom != "X":
			a_to_x += 1
		elif retro_chrom != "X":
			x_to_a += 1
		else:
			x_to_x += 1

a_to_a_pct = a_to_a/total*100
a_to_x_pct = a_to_x/total*100
x_to_a_pct = x_to_a/total*100
x_to_x_pct = x_to_x/total*100

with open(snakemake.output[0], "w") as output_file:
	output_file.write(f"A to A: {a_to_a} ({a_to_a_pct:.2f}%)\n")
	output_file.write(f"A to X: {a_to_x} ({a_to_x_pct:.2f}%)\n")
	output_file.write(f"X to A: {x_to_a} ({x_to_a_pct:.2f}%)\n")
	output_file.write(f"X to X: {x_to_x} ({x_to_x_pct:.2f}%)\n")
