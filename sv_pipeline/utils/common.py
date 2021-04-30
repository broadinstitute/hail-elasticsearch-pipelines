import subprocess
from datetime import datetime

CHROMOSOMES = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21',
    '22', 'X', 'Y', 'M',
]
CHROM_TO_XPOS_OFFSET = {chrom: (1 + i)*int(1e9) for i, chrom in enumerate(CHROMOSOMES)}

GS_SAMPLE_PATH = 'gs://seqr-datasets/v02/GRCh38/RDG_{sample_type}_Broad_Internal/base/projects/{project_guid}/{project_guid}_{file_ext}'


def _get_gs_samples(project_guid, file_ext, sample_type, expected_header):
    """
    Get sample metadata from files in google cloud

    :param project_guid: seqr project identifier
    :param file_ext: extension for the desired sample file
    :param sample_type: sample type (WES/WGS)
    :param expected_header: expected header to validate file
    :return: parsed data from the sample file as a list of lists
    """
    file = GS_SAMPLE_PATH.format(project_guid=project_guid, sample_type=sample_type, file_ext=file_ext)
    process = subprocess.Popen(
        'gsutil cat {}'.format(file), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if process.wait() != 0:
        return None
    header = next(process.stdout).decode('utf-8')
    if header.strip() != expected_header:
        raise Exception('Missing header for sample file, expected "{}" but found {}'.format(
            expected_header, header))
    return [line.decode('utf-8').strip().split('\t') for line in process.stdout]


def get_sample_subset(project_guid, sample_type):
    """
    Get sample id subset for a given project

    :param project_guid: seqr project identifier
    :param sample_type: sample type (WES/WGS)
    :return: set of sample ids
    """
    subset = _get_gs_samples(project_guid, file_ext='ids.txt', sample_type=sample_type, expected_header='s')
    if not subset:
        raise Exception('No sample subset file found')
    return {row[0] for row in subset}


def get_sample_remap(project_guid, sample_type):
    """
    Get an optional remapping for sample ids in the given project

    :param project_guid: seqr project identifier
    :param sample_type: sample type (WES/WGS)
    :return: dictionary mapping VCF sample ids to seqr sample ids, or None if no mapping available
    """
    remap = _get_gs_samples(project_guid, file_ext='remap.tsv', sample_type=sample_type, expected_header='s\tseqr_id')
    if remap:
        remap = {row[0]: row[1] for row in remap}
    return remap


def get_es_index_name(project, meta):
    """
    Get the name for the output ES index

    :param project: seqr project identifier
    :param meta: index metadata
    :return: index name
    """
    return '{project}__structural_variants__{sample_type}__grch{genome_version}__{datestamp}'.format(
        project=project,
        sample_type=meta['sampleType'],
        genome_version=meta['genomeVersion'],
        datestamp=datetime.today().strftime('%Y%m%d'),
    ).lower()
