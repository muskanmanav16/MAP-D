from pathlib import Path

from click.testing import CliRunner
from mapd.cli import get_entity_dict, get_abstract_info, predict_entities, query_db

runner = CliRunner()


class TestCli:
    """Class for testing the CLI commands."""

    def test_get_entity_dict(self):
        result = runner.invoke(get_entity_dict, ['-r'])
        assert result.exit_code == 0

        help_result = runner.invoke(get_entity_dict, ['--help'])
        assert "Enrich the graph with RNA and DNA molecules." in help_result.output

    def test_get_abstract_info(self):
        result = runner.invoke(get_abstract_info, [36316711])
        assert result.exit_code == 0

    def test_predict_entities(self):
        result = runner.invoke(predict_entities, ["Merkel cell polyomavirus (MCPyV) is the only human polyomavirus currently known to cause human cancer."])
        assert result.exit_code == 0

    def test_query_db(self):
        result = runner.invoke(query_db, ['nose', 'test_results.csv', '-s', '2021-01-03', '-e', '2021-09-02'])
        assert result.exit_code == 0


