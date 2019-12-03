from click.testing import CliRunner
from hicreppy.cli import cli

SAMPLE1 = 'data_test/sample_1.cool'
SAMPLE2 = 'data_test/sample_2.cool'

def test_scc():
    """Test exit code and output of scc"""
    runner = CliRunner()
    result = runner.invoke(
        cli, ["scc", SAMPLE1, SAMPLE2]
    )
    assert result.exit_code == 0
    output = float(result.output.split('\n')[-2])
    assert abs(output) < 1

def test_htrain():
    """Test exit code and output of htrain"""
    runner = CliRunner()
    result = runner.invoke(
        cli, ["htrain", SAMPLE1, SAMPLE2]
    )
    assert result.exit_code == 0
    output = int(result.output.split('\n')[-2])
    assert output < 10
    assert output > 0
