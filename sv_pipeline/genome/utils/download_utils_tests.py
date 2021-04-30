import unittest
import mock
import responses

from sv_pipeline.genome.utils.download_utils import download_file

TEST_DIR = 'test/dir'
TEST_TXT_FILE = 'test_file.txt'
TEST_GZ_FILE = 'test_file.gz'
TXT_DATA_URL = 'https://mock_url/test_file.txt'
GZ_DATA_URL = 'https://mock_url/test_file.gz'


class DownloadUtilsTest(unittest.TestCase):

    @responses.activate
    @mock.patch('sv_pipeline.genome.utils.download_utils.logger')
    @mock.patch('sv_pipeline.genome.utils.download_utils.os')
    @mock.patch('sv_pipeline.genome.utils.download_utils.open')
    def test_download_file(self, mock_open, mock_os, mock_logger):
        responses.add(responses.HEAD, GZ_DATA_URL,
                      headers={"Content-Length": "1024"}, status=200)
        responses.add(responses.HEAD, TXT_DATA_URL,
                      headers={"Content-Length": "1024"}, status=200)
        responses.add(responses.GET, TXT_DATA_URL, body='test data\nanother line\n')

        # Test bad url
        with self.assertRaises(ValueError) as ve:
            download_file("bad_url")
        self.assertEqual(str(ve.exception), "Invalid url: bad_url")

        # Test already downloaded
        mock_os.path.isfile.return_value = True
        mock_os.path.getsize.return_value = 1024
        test_file_name = '{}/{}'.format(TEST_DIR, TEST_GZ_FILE)
        mock_os.path.join.return_value = test_file_name
        result = download_file(GZ_DATA_URL)
        self.assertEqual(result, test_file_name)
        mock_logger.info.assert_called_with('Re-using {} previously downloaded from {}'.format(test_file_name, GZ_DATA_URL))

        mock_os.path.isfile.return_value = False
        mock_os.path.getsize.return_value = 0
        test_file_name = "{}/{}".format(TEST_DIR, TEST_TXT_FILE)
        mock_os.path.join.return_value = test_file_name
        mock_os.path.basename.return_value = TEST_TXT_FILE
        mock_logger.reset_mock()
        result = download_file(TXT_DATA_URL, TEST_DIR)
        mock_logger.info.assert_called_with("Downloading {} to {}".format(TXT_DATA_URL, test_file_name))
        self.assertEqual(result, test_file_name)
        mock_os.path.basename.assert_called_with(TXT_DATA_URL)
        mock_os.path.join.assert_called_with(TEST_DIR, TEST_TXT_FILE)
        mock_open.assert_called_with(test_file_name, 'wb')
        handle = mock_open().__enter__()
        handle.writelines.assert_called_once()
