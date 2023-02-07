# from pathlib import Path

from click.testing import CliRunner
from mapd.cli import get_entity_dict, get_abstract_info, predict_entities, query_db

runner = CliRunner()


class TestCli:
    """Class for testing the CLI commands."""

    def test_get_entity_dict(self):
        """Checks if the method to get entities works as expected"""

        result = runner.invoke(get_entity_dict)
        assert result.exit_code == 0

        help_result = runner.invoke(get_entity_dict, ['--help'])
        assert "option to print dict entries row by row" in help_result.output

    def test_get_abstract_info(self):
        """Checks if the method to retrieve abstract information works as expected"""

        # expected_out = {
        #     "pubmed_id": '36298759',
        #     "Title": "Current In Vitro and In Vivo Models to Study MCPyV-Associated MCC.",
        #     "date": "2022-10-07",
        #     "abstract_text": "Merkel cell polyomavirus (MCPyV) is the only human polyomavirus currently known to cause human cancer. MCPyV is believed to be an etiological factor in at least 80% of cases of the rare but aggressive skin malignancy Merkel cell carcinoma (MCC). In these MCPyV+ MCC tumors, clonal integration of the viral genome results in the continued expression of two viral proteins: the viral small T antigen (ST) and a truncated form of the viral large T antigen. The oncogenic potential of MCPyV and the functional properties of the viral T antigens that contribute to neoplasia are becoming increasingly well-characterized with the recent development of model systems that recapitulate the biology of MCPyV+ MCC. In this review, we summarize our understanding of MCPyV and its role in MCC, followed by the current state of both in vitro and in vivo model systems used to study MCPyV and its contribution to carcinogenesis. We also highlight the remaining challenges within the field and the major considerations related to the ongoing development of in vitro and in vivo models of MCPyV+ MCC.",
        #     "entities": {'CANCER': 'skin malignancy Merkel cell carcinoma',
        #       'CELL': 'Merkel cell polyomavirus',
        #       'ORGANISM': 'human polyomavirus'},
        # }
        result = runner.invoke(get_abstract_info, ['36298759'])
        # assert result.output == expected_out
        assert result.exit_code == 0

    def test_predict_entities(self):
        """Checks if the method to predict entities works as expected"""

        result = runner.invoke(predict_entities, ["Merkel cell polyomavirus (MCPyV) is the only human polyomavirus currently known to cause human cancer."])
        assert result.exit_code == 0

    def test_query_db(self):
        """Checks if the method to query the database works as expected"""

        result = runner.invoke(query_db, ['nose', 'test_results.csv', '-s', '2021-01-03', '-e', '2021-09-02'])
        assert result.exit_code == 0
