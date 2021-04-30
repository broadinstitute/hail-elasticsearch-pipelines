import unittest
import mock

from sv_pipeline.genome.utils.mapping_gene_ids import load_gencode, _load_parsed_data_or_download, _parse_gtf_data

DOWNLOAD_PATH = 'test/path'
DOWNLOAD_FILE = 'test/path/gencode.v29.annotation.gtf.gz'
PICKLE_FILE = 'test/path/gencode.v29.annotation.gtf.pickle'
PICKLE_FILE_HANDLE = 'handle'
GTF_DATA = [
    '#description: evidence-based annotation of the human genome, version 31 (Ensembl 97), mapped to GRCh37 with gencode-backmap\n',
    'chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5_2"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2_2"; remap_status "full_contig"; remap_num_mappings 1; remap_target_status "overlap";\n',
    'chr1	HAVANA	gene	621059	622053	.	-	.	gene_id "ENSG00000284662.1_2"; gene_type "protein_coding"; gene_name "OR4F16"; level 2; hgnc_id "HGNC:15079"; havana_gene "OTTHUMG00000002581.3_2"; remap_status "full_contig"; remap_num_mappings 1; remap_target_status "overlap";\n',
    'GL000193.1	HAVANA	gene	77815	78162	.	+	.	gene_id "ENSG00000279783.1_5"; gene_type "processed_pseudogene"; gene_name "AC018692.2"; level 2; havana_gene "OTTHUMG00000189459.1_5"; remap_status "full_contig"; remap_num_mappings 1; remap_target_status "new";\n',
]
GENE_ID_MAPPING = {"DDX11L1": "ENSG00000223972", "OR4F16": "ENSG00000284662", "AC018692.2": "ENSG00000279783"}


class LoadGencodeTestCase(unittest.TestCase):

    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.logger')
    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.pickle')
    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.os.path.isfile')
    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.open')
    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.download_file')
    def _test_load_parsed_data_or_download(self, mock_download, mock_open, mock_isfile, mock_pickle, mock_logger):
        # load from saved pickle
        mock_isfile.return_value = True
        mock_pickle.load.return_value = GENE_ID_MAPPING
        gene_id_mapping, download_file = _load_parsed_data_or_download(29, DOWNLOAD_PATH)
        mock_isfile.assert_called_with(PICKLE_FILE)
        mock_pickle.load.assert_called_with(mock_open.return_value.__enter__.return_value)
        mock_open.assert_called_with(PICKLE_FILE, 'rb')
        self.assertDictEqual(gene_id_mapping, GENE_ID_MAPPING)
        mock_logger.info.assert_called_with('Use the existing pickle file {}.\nIf you want to reload the data, please delete it and re-run the data loading.'.format(PICKLE_FILE))

        # test downloading gtf file
        mock_isfile.return_value = False
        mock_download.return_value = DOWNLOAD_FILE
        mock_logger.reset_mock()
        gene_id_mapping, download_file = _load_parsed_data_or_download(29, DOWNLOAD_PATH)
        mock_isfile.assert_has_calls([mock.call(PICKLE_FILE), mock.call(DOWNLOAD_FILE)])
        mock_download.assert_called_with("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz", to_dir=DOWNLOAD_PATH)
        self.assertEqual(gene_id_mapping, {})
        mock_logger.info.assert_called_with('Downloaded to {}'.format(DOWNLOAD_FILE))
        self.assertEqual(download_file, DOWNLOAD_FILE)

        # test using downloaded file
        mock_isfile.side_effect = [False, True]
        mock_download.reset_mock()
        gene_id_mapping, download_file = _load_parsed_data_or_download(29, DOWNLOAD_PATH)
        mock_isfile.assert_has_calls([mock.call(PICKLE_FILE), mock.call(DOWNLOAD_FILE)])
        mock_logger.info.assert_called_with('Use the existing downloaded file {}. If you want to re-download it, please delete the file and re-run the pipeline.'.format(
            DOWNLOAD_FILE))
        self.assertEqual(download_file, DOWNLOAD_FILE)

    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.pickle')
    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.open')
    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.gzip.open')
    def test_parse_gtf_data(self, mock_gopen, mock_open, mock_pickle):
        # load and parse the gtf file and save to a pickle
        mock_gopen.return_value.__enter__.return_value = GTF_DATA
        gene_id_mapping = _parse_gtf_data(DOWNLOAD_FILE)
        mock_gopen.assert_called_with(DOWNLOAD_FILE, 'rt')
        self.assertEqual(gene_id_mapping, GENE_ID_MAPPING)
        mock_open.assert_called_with(PICKLE_FILE, 'wb')
        mock_pickle.dump.assert_called_with(gene_id_mapping, mock_open.return_value.__enter__.return_value,
                                            protocol=mock_pickle.HIGHEST_PROTOCOL)

        # bad gtf data test
        mock_gopen.return_value.__enter__.return_value = ['bad data']
        with self.assertRaises(ValueError) as ve:
            _ = _parse_gtf_data(DOWNLOAD_FILE)
        self.assertEqual(str(ve.exception), "Unexpected number of fields on line #0: ['bad data']")

    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids.logger')
    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids._load_parsed_data_or_download')
    @mock.patch('sv_pipeline.genome.utils.mapping_gene_ids._parse_gtf_data')
    def test_load_gencode(self, mock_parse_gtf, mock_load, mock_logger):
        # test using saved file
        mock_load.return_value = (GENE_ID_MAPPING, DOWNLOAD_FILE)
        gene_id_mapping = load_gencode(23, download_path=DOWNLOAD_PATH)
        mock_load.assert_called_with(23, DOWNLOAD_PATH)
        mock_logger.info.assert_called_with('Got 3 gene id mapping records')
        mock_parse_gtf.assert_not_called()
        self.assertEqual(gene_id_mapping, GENE_ID_MAPPING)

        # test parsing gtf data
        mock_load.return_value = ({}, DOWNLOAD_FILE)
        mock_parse_gtf.return_value = GENE_ID_MAPPING
        gene_id_mapping = load_gencode(23, download_path=DOWNLOAD_PATH)
        mock_load.assert_called_with(23, DOWNLOAD_PATH)
        mock_parse_gtf.assert_called_with(DOWNLOAD_FILE)
        mock_logger.info.assert_called_with('Got 3 gene id mapping records')
        self.assertEqual(gene_id_mapping, GENE_ID_MAPPING)
