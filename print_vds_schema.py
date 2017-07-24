import argparse
import hail

from pprint import pprint

hc = hail.HailContext()

p = argparse.ArgumentParser()
p.add_argument("input_path", help="input VCF or VDS")
args = p.parse_args()
input_path = args.input_path

print("Input path: %s" % input_path)

if input_path.endswith(".vds"):
    vds = hc.read(input_path)
else:
    vds = hc.import_vcf(input_path, min_partitions=1000, force_bgz=True)

pprint(vds.sample_schema)
pprint(vds.variant_schema)
pprint(vds.genotype_schema)
pprint("%d sample ids" % len(vds.sample_ids))
