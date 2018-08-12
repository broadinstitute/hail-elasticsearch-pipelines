#!/usr/bin/env python

import os

from hail_scripts.v01.utils.add_combined_reference_data import add_combined_reference_data_to_vds
from hail_scripts.v01.utils.add_primate_ai import add_primate_ai_to_vds
from hail_scripts.v01.utils.validate_vds import validate_vds_genome_version_and_sample_type

os.system("pip install elasticsearch")

import argparse
import datetime
import json
import logging
from pprint import pprint
import requests
import time
import sys

from hail_scripts.shared.elasticsearch_utils import ELASTICSEARCH_INDEX, \
    ELASTICSEARCH_UPDATE, ELASTICSEARCH_UPSERT
from hail_scripts.v01.utils.add_hgmd import add_hgmd_to_vds
from hail_scripts.v01.utils.add_eigen import add_eigen_to_vds
from hail_scripts.v01.utils.gcloud_utils import delete_gcloud_file
from hail_scripts.v01.utils.vds_utils import read_in_dataset, compute_minimal_schema, write_vds, run_vep
from hail_scripts.v01.utils.computed_fields import get_expr_for_variant_id, \
    get_expr_for_vep_gene_ids_set, get_expr_for_vep_transcript_ids_set, \
    get_expr_for_vep_consequence_terms_set, \
    get_expr_for_vep_sorted_transcript_consequences_array, \
    get_expr_for_worst_transcript_consequence_annotations_struct, get_expr_for_end_pos, \
    get_expr_for_xpos, get_expr_for_contig, get_expr_for_start_pos, get_expr_for_alt_allele, \
    get_expr_for_ref_allele, get_expr_for_vep_protein_domains_set, get_expr_for_variant_type
from hail_scripts.v01.utils.elasticsearch_utils import DEFAULT_GENOTYPE_FIELDS_TO_EXPORT, \
    ELASTICSEARCH_MAX_SIGNED_SHORT_INT_TYPE, DEFAULT_GENOTYPE_FIELD_TO_ELASTICSEARCH_TYPE_MAP
from hail_scripts.v01.utils.elasticsearch_client import ElasticsearchClient
from hail_scripts.v01.utils.fam_file_utils import MAX_SAMPLES_PER_INDEX, compute_sample_groups_from_fam_file
from hail_scripts.v01.utils.vds_schema_string_utils import convert_vds_schema_string_to_annotate_variants_expr
from hail_scripts.v01.utils.add_1kg_phase3 import add_1kg_phase3_to_vds
from hail_scripts.v01.utils.add_cadd import add_cadd_to_vds
from hail_scripts.v01.utils.add_dbnsfp import add_dbnsfp_to_vds
from hail_scripts.v01.utils.add_clinvar import add_clinvar_to_vds
from hail_scripts.v01.utils.add_exac import add_exac_to_vds
from hail_scripts.v01.utils.add_gnomad import add_gnomad_to_vds
from hail_scripts.v01.utils.add_topmed import add_topmed_to_vds
from hail_scripts.v01.utils.add_mpc import add_mpc_to_vds


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)


p = argparse.ArgumentParser()
p.add_argument("--genome-version", help="Genome build: 37 or 38", choices=["37", "38"], required=True)

p.add_argument("--skip-vep", action="store_true", help="Don't run vep.")
p.add_argument("--skip-annotations", action="store_true", help="Don't add any reference data. Intended for testing.")
p.add_argument('--subset', const="X:31097677-33339441", nargs='?',
               help="All data will first be subsetted to this chrom:start-end range. Intended for testing.")

p.add_argument('--remap-sample-ids', help="Filepath containing 2 tab-separated columns: current sample id and desired sample id")
p.add_argument('--subset-samples', help="Filepath containing ids for samples to keep; if used with --remap-sample-ids, ids are the desired ids (post remapping)")
p.add_argument("--ignore-extra-sample-ids-in-tables", action="store_true")
p.add_argument("--ignore-extra-sample-ids-in-vds", action="store_true")

p.add_argument("--fam-file", help=".fam file used to check VDS sample IDs and assign samples to indices with "
                                  "a max of 'num_samples' per index, but making sure that samples from the same family don't end up in different indices. "
                                  "If used with --remap-sample-ids, contains IDs of samples after remapping")
p.add_argument("--max-samples-per-index", help="Max samples per index", type=int, default=MAX_SAMPLES_PER_INDEX)

p.add_argument('--export-vcf', action="store_true", help="Write out a new VCF file after import")

p.add_argument("--project-guid", help="seqr Project id", required=True)
p.add_argument("--family-id", help="(optional) seqr Family id for datasets (such as Manta SV calls) that are generated per-family")
p.add_argument("--individual-id", help="(optional) seqr Individual id for datasets (such as single-sample Manta SV calls) that are generated per-individual")
p.add_argument("--sample-type", help="sample type (WES, WGS, RNA)", choices=["WES", "WGS", "RNA"], required=True)
p.add_argument("--dataset-type", help="what pipeline was used to generate the data", choices=["GATK_VARIANTS", "MANTA_SVS", "JULIA_SVS"], required=True)

p.add_argument("--index", help="(optional) elasticsearch index name. If not specified, the index name will be computed based on project_guid, family_id, sample_type and dataset_type.")

p.add_argument("--host", help="Elastisearch IP address", default="10.56.10.4")
p.add_argument("--port", help="Elastisearch port", default="9200")
p.add_argument("--num-shards", help="Number of index shards", type=int, default=12)

p.add_argument("--vep-block-size", help="Block size to use for VEP", default=200, type=int)
p.add_argument("--es-block-size", help="Block size to use when exporting to elasticsearch", default=1000, type=int)


p.add_argument("--exclude-dbnsfp", action="store_true", help="Don't add annotations from dbnsfp. Intended for testing.")
p.add_argument("--exclude-1kg", action="store_true", help="Don't add 1kg AFs. Intended for testing.")
p.add_argument("--exclude-omim", action="store_true", help="Don't add OMIM mim id column. Intended for testing.")
p.add_argument("--exclude-gene-constraint", action="store_true", help="Don't add gene constraint columns. Intended for testing.")
p.add_argument("--exclude-eigen", action="store_true", help="Don't add Eigen scores. Intended for testing.")
p.add_argument("--exclude-cadd", action="store_true", help="Don't add CADD scores (they take a really long time to load). Intended for testing.")
p.add_argument("--exclude-gnomad", action="store_true", help="Don't add gnomAD exome or genome fields. Intended for testing.")
p.add_argument("--exclude-exac", action="store_true", help="Don't add ExAC fields. Intended for testing.")
p.add_argument("--exclude-topmed", action="store_true", help="Don't add TopMed AFs. Intended for testing.")
p.add_argument("--exclude-clinvar", action="store_true", help="Don't add clinvar fields. Intended for testing.")
p.add_argument("--exclude-hgmd", action="store_true", help="Don't add HGMD fields. Intended for testing.")
p.add_argument("--exclude-mpc", action="store_true", help="Don't add MPC fields. Intended for testing.")
p.add_argument("--exclude-primate-ai", action="store_true", help="Don't add PrimateAI fields. Intended for testing.")
p.add_argument("--exclude-gnomad-coverage", action="store_true", help="Don't add gnomAD exome and genome coverage. Intended for testing.")
p.add_argument("--exclude-vcf-info-field", action="store_true", help="Don't add any fields from the VCF info field. Intended for testing.")

p.add_argument("--dont-update-operations-log", action="store_true", help="Don't save metadata about this export in the operations log.")
p.add_argument("--create-snapshot", action="store_true", help="Create an elasticsearch snapshot in a google bucket after indexing is complete.")
p.add_argument("--dont-delete-intermediate-vds-files", action="store_true", help="Keep intermediate VDS files to allow restarting the pipeline from the middle using --start-with-step")

p.add_argument("--start-with-step", help="Which pipeline step to start with.", type=int, default=0, choices=[0, 1, 2, 3, 4])
p.add_argument("--start-with-sample-group", help="If the callset contains more samples than the limit specified by --max-samples-per-index, "
                                                 "it will be loaded into multiple separate indices. Setting this command-line arg to a value > 0 causes the pipeline to start from sample "
                                                 "group other than the 1st one. This is useful for restarting a failed pipeline from exactly where it left off.", type=int, default=0)

p.add_argument("--username", help="(optional) user running this pipeline. This is the local username and it must be passed in because the script can't look it up when it runs on dataproc.")
p.add_argument("--directory", help="(optional) current directory. This is the local directory and it must be passed in because the script can't look it up when it runs on dataproc.")

p.add_argument("--output-vds", help="(optional) Output vds filename prefix (eg. test-vds)")

p.add_argument("input_vds", help="input VDS")

args = p.parse_args()

if args.dataset_type == "GATK_VARIANTS":
    variant_type_string = "variants"
elif args.dataset_type in ["MANTA_SVS", "JULIA_SVS"]:
    variant_type_string = "sv"
else:
    raise ValueError("Unexpected args.dataset_type == " + str(args.dataset_type))


# generate the index name as:  <project>_<WGS_WES>_<family?>_<VARIANTS or SVs>_<YYYYMMDD>_<batch>
if args.index:
    index_name = args.index.lower()
else:
    index_name = "%s%s%s__%s__grch%s__%s__%s" % (
        args.project_guid,
        "__"+args.family_id if args.family_id else "",  # optional family id
        "__"+args.individual_id if args.individual_id else "",  # optional individual id
        args.sample_type,
        args.genome_version,
        variant_type_string,
        datetime.datetime.now().strftime("%Y%m%d"),
    )

    index_name = index_name.lower()  # elasticsearch requires index names to be all lower-case

logger.info("Index name: %s" % (index_name,))


def export_to_elasticsearch(
    host,
    port,
    vds,
    index_name,
    args,
    operation=ELASTICSEARCH_INDEX,
    delete_index_before_exporting=False,
    export_genotypes=True,
    disable_doc_values_for_fields=(),
    disable_index_for_fields=(),
    export_snapshot_to_google_bucket=False,
    update_operations_log=False,
    start_with_sample_group=0,
):
    """Utility method for exporting the given vds to an elasticsearch index.
    """

    logger.info("Input: " + input_path)

    index_type = "variant"

    if export_genotypes:
        if args.dataset_type == "GATK_VARIANTS":
            genotype_fields_to_export = DEFAULT_GENOTYPE_FIELDS_TO_EXPORT
            genotype_field_to_elasticsearch_type_map = DEFAULT_GENOTYPE_FIELD_TO_ELASTICSEARCH_TYPE_MAP
        elif args.dataset_type in ["MANTA_SVS", "JULIA_SVS"]:
            genotype_fields_to_export = [
                'num_alt = if(g.GT.isCalled()) g.GT.nNonRefAlleles() else -1',
                #'genotype_filter = g.FT',
                #'gq = g.GQ',
                'dp = if(g.GT.isCalled()) [g.PR.sum + g.SR.sum, '+ELASTICSEARCH_MAX_SIGNED_SHORT_INT_TYPE+'].min() else NA:Int',
                'ab = let total=g.PR.sum + g.SR.sum in if(g.GT.isCalled() && total != 0) ((g.PR[1] + g.SR[1]) / total).toFloat else NA:Float',
                'ab_PR = let total=g.PR.sum in if(g.GT.isCalled() && total != 0) (g.PR[1] / total).toFloat else NA:Float',
                'ab_SR = let total=g.SR.sum in if(g.GT.isCalled() && total != 0) (g.SR[1] / total).toFloat else NA:Float',
                'dp_PR = if(g.GT.isCalled()) [g.PR.sum,'+ELASTICSEARCH_MAX_SIGNED_SHORT_INT_TYPE+'].min() else NA:Int',
                'dp_SR = if(g.GT.isCalled()) [g.SR.sum,'+ELASTICSEARCH_MAX_SIGNED_SHORT_INT_TYPE+'].min() else NA:Int',
            ]

            genotype_field_to_elasticsearch_type_map = {
                ".*_num_alt": {"type": "byte", "doc_values": "false"},
                #".*_genotype_filter": {"type": "keyword", "doc_values": "false"},
                #".*_gq": {"type": "short", "doc_values": "false"},
                ".*_dp": {"type": "short", "doc_values": "false"},
                ".*_ab": {"type": "half_float", "doc_values": "false"},
                ".*_ab_PR": {"type": "half_float", "doc_values": "false"},
                ".*_ab_SR": {"type": "half_float", "doc_values": "false"},
                ".*_dp_PR": {"type": "short", "doc_values": "false"},
                ".*_dp_SR": {"type": "short", "doc_values": "false"},
            }
        else:
            raise ValueError("Unexpected args.dataset_type: %s" % args.dataset_type)
    else:
        genotype_fields_to_export = []
        genotype_field_to_elasticsearch_type_map = {}

    vds = vds.persist()

    client = ElasticsearchClient(host, port)
    for i, sample_group in enumerate(sample_groups):

        if i < start_with_sample_group:
            continue

        #if delete_index_before_exporting and i < 4:
        #    continue

        if len(sample_groups) > 1:
            vds_sample_subset = vds.filter_samples_list(sample_group, keep=True)
            current_index_name = "%s_%s" % (index_name, i)
        else:
            vds_sample_subset = vds
            current_index_name = index_name

        logger.info("==> exporting %s samples into %s" % (len(sample_group), current_index_name))
        logger.info("Samples: %s .. %s" % (", ".join(sample_group[:3]), ", ".join(sample_group[-3:])))

        logger.info("==> export to elasticsearch")
        pprint(vds.variant_schema)

        timestamp1 = time.time()

        client.export_vds_to_elasticsearch(
            vds_sample_subset,
            genotype_fields_to_export=genotype_fields_to_export,
            genotype_field_to_elasticsearch_type_map=genotype_field_to_elasticsearch_type_map,
            index_name=current_index_name,
            index_type_name=index_type,
            block_size=args.es_block_size,
            num_shards=args.num_shards,
            delete_index_before_exporting=delete_index_before_exporting,
            elasticsearch_write_operation=operation,
            elasticsearch_mapping_id="docId",
            disable_doc_values_for_fields=disable_doc_values_for_fields,
            disable_index_for_fields=disable_index_for_fields,
            is_split_vds=True,
            verbose=True,
        )

        timestamp2 = time.time()
        logger.info("==> finished exporting - time: %s seconds" % (timestamp2 - timestamp1))

    if export_snapshot_to_google_bucket:
        logger.info("==> export snapshot to google bucket")
        client.create_elasticsearch_snapshot(
            index_name=index_name + "*",
            bucket="seqr-database-backups",
            base_path="elasticsearch/snapshots",
            snapshot_repo="callsets")

    if update_operations_log:
        logger.info("==> update operations log")
        client.save_index_operation_metadata(
            args.input_vds,
            index_name,
            args.genome_version,
            fam_file=args.fam_file,
            remap_sample_ids=args.remap_sample_ids,
            subset_samples=args.subset_samples,
            skip_vep=args.skip_vep,
            project_id=args.project_guid,
            dataset_type=args.dataset_type,
            sample_type=args.sample_type,
            command=" ".join(sys.argv),
            directory=args.directory,
            username=args.username,
            operation="create_index",
            status="success",
        )


input_path = str(args.input_vds).rstrip("/")
if not (input_path.endswith(".vds") or input_path.endswith(".vcf") or input_path.endswith(".vcf.gz") or input_path.endswith(".vcf.bgz")):
    p.error("Input must be a .vds or .vcf.gz")

input_path_prefix = input_path.replace(".vds", "")

filter_interval = "1-MT"
if args.subset:
    filter_interval = args.subset


logger.info("\n==> create HailContext")

import hail  # import hail here so that you can run this script with --help even if hail isn't installed locally.
hc = hail.HailContext(log="/hail.log")

logger.info("Reading in dataset...")
vds = read_in_dataset(hc, input_path, dataset_type=args.dataset_type, filter_interval=filter_interval)

validate_vds_genome_version_and_sample_type(hc, vds, args.genome_version, args.sample_type)


output_vds_hash = ""

# NOTE: if sample IDs are remapped first thing, then the fam file should contain the desired (not original IDs)
if args.remap_sample_ids:
    logger.info("Remapping sample ids...")
    id_map = hc.import_table(args.remap_sample_ids, no_header=True)
    mapping = dict(zip(id_map.query('f0.collect()'), id_map.query('f1.collect()')))
    # check that ids being remapped exist in VDS
    samples_in_table = set(mapping.keys())
    samples_in_vds = set(vds.sample_ids)
    matched = samples_in_table.intersection(samples_in_vds)
    if len(matched) < len(samples_in_table):
        warning_message = ("Only {0} out of {1} remapping-table IDs matched IDs in the variant callset.\n"
            "Remapping-table IDs that aren't in the VDS: {2}\n"
            "All VDS IDs: {3}").format(
            len(matched), len(samples_in_table), list(samples_in_table.difference(samples_in_vds)), samples_in_vds)
        if not args.ignore_extra_sample_ids_in_tables:
            raise ValueError(warning_message)
        logger.warning(warning_message)
    vds = vds.rename_samples(mapping)
    logger.info('Remapped {} sample ids...'.format(len(matched)))


# subset samples as desired
if args.subset_samples:
    logger.info("Subsetting to specified samples...")
    keep_samples = hc.import_table(args.subset_samples, no_header=True).key_by('f0')
    # check that all subset samples exist in VDS
    samples_in_table = set(keep_samples.query('f0.collect()'))
    samples_in_vds = set(vds.sample_ids)
    matched = samples_in_table.intersection(samples_in_vds)
    if len(matched) < len(samples_in_table):
        warning_message = ("Only {0} out of {1} subsetting-table IDs matched IDs in the variant callset.\n" \
            "Dropping {2} IDs that aren't in the VDS: {3}\n"
            "All VDS IDs: {4}").format(
            len(matched), len(samples_in_table), len(samples_in_table) - len(matched), list(samples_in_table.difference(samples_in_vds)), samples_in_vds)
        if not args.ignore_extra_sample_ids_in_tables:
            raise ValueError(warning_message)
        logger.warning(warning_message)
    original_sample_count = vds.num_samples
    vds = vds.filter_samples_table(keep_samples, keep=True).variant_qc().filter_variants_expr('va.qc.AC > 0')
    new_sample_count = vds.num_samples
    logger.info('Kept {0} out of {1} samples in vds'.format(new_sample_count, original_sample_count))

    output_vds_hash = "_%s_samples__%020d" % (len(matched), abs(hash(",".join(sorted(list(matched))))))

    logger.info("Finished subsetting samples.")
    logger.info("Callset stats after subsetting:")
    summary = vds.summarize()
    pprint(summary)

# compute sample groups
if len(vds.sample_ids) > args.max_samples_per_index:
    if not args.fam_file:
        p.exit("--fam-file must be specified for callsets larger than %s samples. This callset has %s samples." % (args.max_samples_per_index, len(vds.sample_ids)))
    else:
        sample_groups = compute_sample_groups_from_fam_file(
            args.fam_file,
            vds.sample_ids,
            args.max_samples_per_index,
            args.ignore_extra_sample_ids_in_vds,
            args.ignore_extra_sample_ids_in_tables,
        )
else:
    sample_groups = [vds.sample_ids]

if args.output_vds:
    output_vds_prefix = os.path.join(os.path.dirname(input_path), args.output_vds.replace(".vds", ""))
else:
    output_vds_prefix = input_path.replace(".vcf", "").replace(".vds", "").replace(".bgz", "").replace(".gz", "") + output_vds_hash

step0_output_vds = output_vds_prefix + (".vep.vds" if ".vep" not in output_vds_prefix and not args.skip_vep else ".vds")
step1_output_vds = output_vds_prefix + ".vep_and_computed_annotations.vds"
step3_output_vds = output_vds_prefix + ".vep_and_all_annotations.vds"

# Store step0_output_vds as the cached version of the dataset in google buckets, and also set it as the global.sourceFilePath
# because
# 1) vep is the most time-consuming step (other than exporting to elasticsearch), so it makes sense to cache results
# 2) at this stage, all subsetting and remapping has already been applied, so the samples in the dataset are only the ones exported to elasticsearch
# 3) annotations may be updated / added more often than vep versions.

vds = vds.annotate_global_expr('global.sourceFilePath = "{}"'.format(step0_output_vds))

# run vep
if args.start_with_step == 0:
    if not args.skip_vep:
        logger.info("\n\n=============================== pipeline - step 0 ===============================")
        logger.info("Read in data, run vep, write data to VDS")
        vds = run_vep(vds, genome_version=args.genome_version, block_size=args.vep_block_size)

        vds = vds.annotate_global_expr('global.gencodeVersion = "{}"'.format("19" if args.genome_version == "37" else "25"))

    if step0_output_vds != input_path:
        write_vds(vds, step0_output_vds)

    # write out new vcf (after sample id remapping and subsetting if requested above)
    if args.export_vcf:
        logger.info("Writing out to VCF...")
        if not args.skip_vep:
            vds.export_vcf(output_vds_prefix + ".vep.vcf.bgz", overwrite=True)
        else:
            vds.export_vcf(output_vds_prefix + ".vcf.bgz", overwrite=True)

hc.stop()


if args.start_with_step <= 1:
    logger.info("\n\n=============================== pipeline - step 1 ===============================")
    logger.info("Read in data, compute various derived fields, write data to VDS")

    logger.info("\n==> re-create HailContext")
    hc = hail.HailContext(log="/hail.log")

    vds = read_in_dataset(hc, step0_output_vds, dataset_type=args.dataset_type, filter_interval=filter_interval, skip_summary=True)

    # add computed annotations
    logger.info("\n==> adding computed annotations")
    parallel_computed_annotation_exprs = [
        "va.docId = %s" % get_expr_for_variant_id(512),
        "va.variantId = %s" % get_expr_for_variant_id(),

        "va.variantType= %s" % get_expr_for_variant_type(),
        "va.contig = %s" % get_expr_for_contig(),
        "va.pos = %s" % get_expr_for_start_pos(),
        "va.start = %s" % get_expr_for_start_pos(),
        "va.end = %s" % get_expr_for_end_pos(),
        "va.ref = %s" % get_expr_for_ref_allele(),
        "va.alt = %s" % get_expr_for_alt_allele(),
        "va.xpos = %s" % get_expr_for_xpos(pos_field="start"),
        "va.xstart = %s" % get_expr_for_xpos(pos_field="start"),

        "va.transcriptIds = %s" % get_expr_for_vep_transcript_ids_set(vep_transcript_consequences_root="va.vep.transcript_consequences"),
        "va.domains = %s" % get_expr_for_vep_protein_domains_set(vep_transcript_consequences_root="va.vep.transcript_consequences"),
        "va.transcriptConsequenceTerms = %s" % get_expr_for_vep_consequence_terms_set(vep_transcript_consequences_root="va.vep.transcript_consequences"),
        "va.sortedTranscriptConsequences = %s" % get_expr_for_vep_sorted_transcript_consequences_array(vep_root="va.vep"),
    ]

    serial_computed_annotation_exprs = [
        "va.xstop = %s" % get_expr_for_xpos(field_prefix="va.", pos_field="end"),
        "va.mainTranscript = %s" % get_expr_for_worst_transcript_consequence_annotations_struct("va.sortedTranscriptConsequences"),
        "va.geneIds = %s" % get_expr_for_vep_gene_ids_set(vep_transcript_consequences_root="va.sortedTranscriptConsequences", exclude_upstream_downstream_genes=True),
        "va.codingGeneIds = %s" % get_expr_for_vep_gene_ids_set(vep_transcript_consequences_root="va.sortedTranscriptConsequences", only_coding_genes=True, exclude_upstream_downstream_genes=True),
        "va.sortedTranscriptConsequences = json(va.sortedTranscriptConsequences)",
    ]

    vds = vds.annotate_variants_expr(parallel_computed_annotation_exprs)

    for expr in serial_computed_annotation_exprs:
        vds = vds.annotate_variants_expr(expr)

    pprint(vds.variant_schema)

    # apply schema to dataset
    INPUT_SCHEMA  = {}
    if args.dataset_type == "GATK_VARIANTS":
        INPUT_SCHEMA["top_level_fields"] = """
            docId: String,
            variantId: String,
            originalAltAlleles: Set[String],

            contig: String,
            start: Int,
            pos: Int,
            end: Int,
            ref: String,
            alt: String,

            xpos: Long,
            xstart: Long,
            xstop: Long,

            rsid: String,
            --- qual: Double,
            filters: Set[String],
            wasSplit: Boolean,
            aIndex: Int,

            geneIds: Set[String],
            transcriptIds: Set[String],
            codingGeneIds: Set[String],
            domains: Set[String],
            transcriptConsequenceTerms: Set[String],
            sortedTranscriptConsequences: String,
            mainTranscript: Struct,
        """

        INPUT_SCHEMA["info_fields"] = """
            AC: Array[Int],
            AF: Array[Double],
            AN: Int,
            --- BaseQRankSum: Double,
            --- ClippingRankSum: Double,
            --- DP: Int,
            --- FS: Double,
            InbreedingCoeff: Double,
            MQ: Double,
            --- MQRankSum: Double,
            QD: Double,
            --- ReadPosRankSum: Double,
            --- VQSLOD: Double,
            --- culprit: String,
        """
    elif args.dataset_type in ["MANTA_SVS", "JULIA_SVS"]:
        INPUT_SCHEMA["top_level_fields"] = """
            docId: String,
            variantId: String,

            contig: String,
            start: Int,
            pos: Int,
            end: Int,
            ref: String,
            alt: String,

            xpos: Long,
            xstart: Long,
            xstop: Long,

            rsid: String,
            --- qual: Double,
            filters: Set[String],

            geneIds: Set[String],
            transcriptIds: Set[String],
            codingGeneIds: Set[String],
            domains: Set[String],
            transcriptConsequenceTerms: Set[String],
            sortedTranscriptConsequences: String,
            mainTranscript: Struct,
        """

        # END=100371979;SVTYPE=DEL;SVLEN=-70;CIGAR=1M70D	GT:FT:GQ:PL:PR:SR
        if args.dataset_type == "MANTA_SVS":
            INPUT_SCHEMA["info_fields"] = """
                IMPRECISE: Boolean,
                SVTYPE: String,
                SVLEN: Int,
                END: Int,
            """
        else:
            INPUT_SCHEMA["info_fields"] = """
                IMPRECISE: Boolean,
                SVTYPE: String,
                SVLEN: Int,
                END: Int,
                OCC: Int,
                FRQ: Double,
            """
    else:
        raise ValueError("Unexpected dataset_type: %s" % args.dataset_type)

    if args.exclude_vcf_info_field:
        INPUT_SCHEMA["info_fields"] = ""

    expr = convert_vds_schema_string_to_annotate_variants_expr(root="va.clean", **INPUT_SCHEMA)

    vds = vds.annotate_variants_expr(expr=expr)
    vds = vds.annotate_variants_expr("va = va.clean")

    write_vds(vds, step1_output_vds)

    hc.stop()


if args.start_with_step <= 2:
    logger.info("\n\n=============================== pipeline - step 2 ===============================")
    logger.info("Read in data, add more reference datasets, export to elasticsearch")

    logger.info("\n==> create HailContext")
    hc = hail.HailContext(log="/hail.log")

    vds = read_in_dataset(hc, step1_output_vds, dataset_type=args.dataset_type, filter_interval=filter_interval, skip_summary=True)

    #if args.dataset_type == "GATK_VARIANTS":

        #if not args.skip_annotations and not args.exclude_omim:
        #    logger.info("\n==> Add omim info")
        #    vds = add_omim_to_vds(hc, vds, root="va.omim", vds_key='va.mainTranscript.gene_id')

        #if not args.skip_annotations and not args.exclude_gene_constraint:
        #    logger.info("\n==> Add gene constraint")
        #    vds = add_gene_constraint_to_vds(hc, vds)

    export_to_elasticsearch(
        args.host,
        args.port,
        vds,
        index_name,
        args,
        operation=ELASTICSEARCH_UPSERT,
        delete_index_before_exporting=True,
        export_genotypes=True,
        disable_doc_values_for_fields=("sortedTranscriptConsequences", ),
        disable_index_for_fields=("sortedTranscriptConsequences", ),
        export_snapshot_to_google_bucket=False,
        update_operations_log=False,
        start_with_sample_group=args.start_with_sample_group if args.start_with_step == 0 else 0,
    )

    hc.stop()

if args.start_with_step <= 3:
    logger.info("\n\n=============================== pipeline - step 3 ===============================")
    logger.info("Read in data, add more reference datasets, write data to VDS")

    logger.info("\n==> create HailContext")
    hc = hail.HailContext(log="/hail.log")

    vds = read_in_dataset(hc, step1_output_vds, dataset_type=args.dataset_type, filter_interval=filter_interval, skip_summary=True)
    vds = compute_minimal_schema(vds, args.dataset_type)

    if args.dataset_type == "GATK_VARIANTS":
        #if not args.skip_annotations and not args.exclude_gnomad_coverage:
        #    logger.info("\n==> Add gnomad coverage")
        #    vds = add_gnomad_exome_coverage_to_vds(hc, vds, args.genome_version, root="va.gnomad_exome_coverage")
        #    vds = add_gnomad_genome_coverage_to_vds(hc, vds, args.genome_version, root="va.gnomad_genome_coverage")
        if not (args.exclude_dbnsfp or args.exclude_cadd or args.exclude_1kg or args.exclude_exac or
                args.exclude_topmed or args.exclude_mpc or args.exclude_gnomad or args.exclude_eigen):
            # annotate with the combined reference data file which was generated using
            # ../download_and_create_reference_datasets/create_all_variant_level_reference_data_vds.py
            # and contains all these annotations in one .vds
            logger.info("\n==> add combined variant-level reference data")
            vds = add_combined_reference_data_to_vds(hc, vds, args.genome_version, subset=filter_interval)

        else:
            # annotate with each reference data file - one-by-one
            if not args.skip_annotations and not args.exclude_dbnsfp:
                logger.info("\n==> add dbnsfp")
                vds = add_dbnsfp_to_vds(hc, vds, args.genome_version, root="va.dbnsfp", subset=filter_interval)

            if not args.skip_annotations and not args.exclude_cadd:
                logger.info("\n==> add cadd")
                vds = add_cadd_to_vds(hc, vds, args.genome_version, root="va.cadd", subset=filter_interval)

            if not args.skip_annotations and not args.exclude_1kg:
                logger.info("\n==> add 1kg")
                vds = add_1kg_phase3_to_vds(hc, vds, args.genome_version, root="va.g1k", subset=filter_interval)

            if not args.skip_annotations and not args.exclude_exac:
                logger.info("\n==> add exac")
                vds = add_exac_to_vds(hc, vds, args.genome_version, root="va.exac", subset=filter_interval)

            if not args.skip_annotations and not args.exclude_topmed:
                logger.info("\n==> add topmed")
                vds = add_topmed_to_vds(hc, vds, args.genome_version, root="va.topmed", subset=filter_interval)

            if not args.skip_annotations and not args.exclude_mpc:
                logger.info("\n==> add mpc")
                vds = add_mpc_to_vds(hc, vds, args.genome_version, root="va.mpc", subset=filter_interval)

            if not args.skip_annotations and not args.exclude_gnomad:
                logger.info("\n==> add gnomad exomes")
                vds = add_gnomad_to_vds(hc, vds, args.genome_version, exomes_or_genomes="exomes", root="va.gnomad_exomes", subset=filter_interval)

            if not args.skip_annotations and not args.exclude_gnomad:
                logger.info("\n==> add gnomad genomes")
                vds = add_gnomad_to_vds(hc, vds, args.genome_version, exomes_or_genomes="genomes", root="va.gnomad_genomes", subset=filter_interval)

            if not args.skip_annotations and not args.exclude_eigen:
                logger.info("\n==> add eigen")
                vds = add_eigen_to_vds(hc, vds, args.genome_version, root="va.eigen", subset=filter_interval)

            if not args.exclude_primate_ai:
                logger.info("\n==> add primate_ai")
                vds = add_primate_ai_to_vds(hc, vds, args.genome_version, root="va.primate_ai", subset=filter_interval)

    if not args.skip_annotations and not args.exclude_clinvar:
        logger.info("\n==> add clinvar")
        vds = add_clinvar_to_vds(hc, vds, args.genome_version, root="va.clinvar", subset=filter_interval)

    if not args.skip_annotations and not args.exclude_hgmd:
        logger.info("\n==> add hgmd")
        vds = add_hgmd_to_vds(hc, vds, args.genome_version, root="va.hgmd", subset=filter_interval)

    write_vds(vds, step3_output_vds)

    hc.stop()


if args.start_with_step <= 4:
    logger.info("\n\n=============================== pipeline - step 4 ===============================")
    logger.info("Read in data, export data to elasticsearch")

    logger.info("\n==> create HailContext")
    hc = hail.HailContext(log="/hail.log")

    vds = read_in_dataset(hc, step3_output_vds, dataset_type=args.dataset_type, filter_interval=filter_interval, skip_summary=True)

    export_to_elasticsearch(
        args.host,
        args.port,
        vds,
        index_name,
        args,
        operation=ELASTICSEARCH_UPDATE,
        delete_index_before_exporting=False,
        export_genotypes=False,
        disable_doc_values_for_fields=(),
        disable_index_for_fields=(),
        export_snapshot_to_google_bucket=args.create_snapshot,
        update_operations_log=not args.dont_update_operations_log,
        start_with_sample_group=args.start_with_sample_group if args.start_with_step == 1 else 0,
    )

    hc.stop()


if not args.dont_delete_intermediate_vds_files:
    #delete_gcloud_file(step0_output_vds) -- don't delete the vep-annotated vds
    delete_gcloud_file(step1_output_vds)
    delete_gcloud_file(step3_output_vds)
