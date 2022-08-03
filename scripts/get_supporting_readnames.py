from collections import defaultdict

read_threshold = snakemake.config["read_threshold"]

transcript_junction_to_reads = defaultdict(set)
with open(snakemake.input[0], "r") as input_spannedjunctions_file:
  for line in input_spannedjunctions_file:
    geneid, transcript, junction, readname = line.strip().split()
    transcript_junction_to_reads[(transcript, junction)].add(readname)

supporting_reads = set()

for (transcript, junction) in transcript_junction_to_reads:
  reads = transcript_junction_to_reads[(transcript, junction)]
  num_reads = len(reads)
  if num_reads >= read_threshold:
    supporting_reads.update(reads)

with open(snakemake.output[0], "w") as output_readnames_file:
  for read in supporting_reads:
    output_readnames_file.write(f"{read}\n")
