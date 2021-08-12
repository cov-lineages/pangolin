from pathlib import Path
from pangolin import command

TEST_DIR = Path(__file__).parent


def test_cmd_line(tmp_path):
    query_file = TEST_DIR / 'test-data' / 'sequence1.fasta'
    output_file = tmp_path / 'out.csv'
    args = ['--outfile', str(output_file), str(query_file)]
    command.main(sysargs=args)
    results = open(output_file).read()
    assert 'Delta (B.1.617.2-like)' in results


def test_cmd_line_partial_datadir(tmp_path):
    # test with a datadir with just the `constellations` data
    test_data_dir = TEST_DIR / 'test-data'
    query_file = test_data_dir / 'sequence1.fasta'
    datadir = test_data_dir / 'datadir1'

    output_file = tmp_path / 'out.csv'
    args = ['--outfile', str(output_file),
            '--datadir', str(datadir), str(query_file)]
    command.main(sysargs=args)
    results = open(output_file).read()
    assert 'Delta (B.1.617.2-like-with-bells)' in results
