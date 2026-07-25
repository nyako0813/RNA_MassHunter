from pathlib import Path
import re

import yaml

import rna_masshunter.sciex_cross_layer_bundle_cli as cli
from rna_masshunter.sciex_cross_layer_bundle_runner import XL_SHEET_NAMES


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
WORKFLOW = ROOT / "docs" / "sciex_cross_layer_bundle_workflow.md"
CI = ROOT / ".github" / "workflows" / "tests.yml"
LOCK = ROOT / "requirements-lock.txt"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_readme_cli_options_match_parser():
    readme = _read(README)
    parser_options = {
        option
        for action in cli.build_parser()._actions
        for option in action.option_strings
        if option.startswith("--") and option != "--help"
    }
    expected = {
        "--config", "--full-bundle", "--t1-bundle", "--p1ap-ms1-bundle",
        "--p1ap-ms2-bundle", "--aggregate-json", "--excel-output",
        "--summary-json", "--dry-run", "--debug",
    }
    assert parser_options == expected
    assert all(option in readme for option in expected)


def test_documented_config_example_loads_with_supported_keys(tmp_path):
    readme = _read(README)
    match = re.search(
        r"<!-- sciex-cross-layer-config-example:start -->\s*```yaml\n(.*?)```\s*"
        r"<!-- sciex-cross-layer-config-example:end -->",
        readme,
        flags=re.DOTALL,
    )
    assert match is not None
    payload = yaml.safe_load(match.group(1))
    assert set(payload) == {"cross_layer_bundle_runner"}
    section = payload["cross_layer_bundle_runner"]
    assert set(section) == cli._CONFIG_KEYS
    assert set(section["bundles"]) == cli._CONFIG_BUNDLE_KEYS
    assert set(section["output"]) == cli._CONFIG_OUTPUT_KEYS
    config_path = tmp_path / "example.yaml"
    config_path.write_text(match.group(1), encoding="utf-8")
    loaded = cli.load_cross_layer_bundle_cli_config(config_path)
    assert loaded.enabled is True
    assert loaded.overwrite is False
    assert all(path is not None and path.is_absolute() for path in loaded.bundles.values())


def test_documented_sheet_names_match_runner_contract():
    readme = _read(README)
    workflow = _read(WORKFLOW)
    assert len(XL_SHEET_NAMES) == 6
    assert all(len(name) <= 31 for name in XL_SHEET_NAMES)
    for name in XL_SHEET_NAMES:
        assert f"`{name}`" in readme
        assert f"`{name}`" in workflow


def test_documented_exit_codes_match_cli_constants():
    readme = _read(README)
    codes = {
        cli.EXIT_SUCCESS,
        cli.EXIT_ARGUMENT_CONFIG,
        cli.EXIT_BUNDLE_VALIDATION,
        cli.EXIT_COMPATIBILITY,
        cli.EXIT_OUTPUT,
        cli.EXIT_AGGREGATION,
    }
    documented = {
        int(value)
        for value in re.findall(r"^\| (\d+) \| .* \|$", readme, flags=re.MULTILINE)
    }
    assert documented == codes


def test_local_documentation_links_resolve():
    for source in (README, WORKFLOW):
        for target in re.findall(r"\[[^]]+\]\(([^)]+)\)", _read(source)):
            if "://" in target or target.startswith("#"):
                continue
            path_text = target.split("#", 1)[0]
            assert (source.parent / path_text).resolve().is_file(), (source, target)


def test_ci_and_dependency_snapshot_are_reproducible_and_synthetic_only():
    ci = _read(CI)
    assert 'python-version: "3.12"' in ci
    assert "pip install -r requirements-lock.txt" in ci
    assert "PYTHONPATH: ." in ci
    assert "pytest -q" in ci
    assert "data/reference" not in ci
    requirements = [
        line for line in _read(LOCK).splitlines()
        if line and not line.startswith("#")
    ]
    assert requirements
    assert all(re.fullmatch(r"[A-Za-z0-9_.-]+==[^=\s]+", line) for line in requirements)
    for direct in ("PyYAML", "pandas", "openpyxl", "numpy", "scipy", "pyteomics", "tqdm", "lxml", "psims", "pytest"):
        assert any(line.lower().startswith(direct.lower() + "==") for line in requirements)


def test_new_text_artifacts_have_no_trailing_whitespace():
    for path in (README, WORKFLOW, CI, LOCK, Path(__file__)):
        bad_lines = [
            number for number, line in enumerate(_read(path).splitlines(), start=1)
            if line != line.rstrip()
        ]
        assert not bad_lines, f"{path}: trailing whitespace on {bad_lines}"
