from pathlib import Path
import pytest
import yaml
from main import parse_args, resolve_config_path
from rna_masshunter.config import load_config, validate_config


def minimal(enzyme="RNase_T1"):
    return {"input":{"mzml_path":"sample.mzML"},"sequence":{"name":"test","sequence":"ACGU","anticodon":"ACG","wobble_position":2},"digestion":{"enabled":True,"enzyme":enzyme},"reconstruction":{"enabled":False}}


def test_default_config_path_is_repository_config():
    root=Path('/repo')
    assert resolve_config_path(root,parse_args([]).config)==root/'config.yaml'


def test_relative_config_path_uses_current_working_directory(tmp_path,monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert resolve_config_path(Path('/repo'),'alternate.yaml')==tmp_path/'alternate.yaml'


def test_absolute_config_path_is_preserved(tmp_path):
    path=tmp_path/'alternate.yaml'
    assert resolve_config_path(Path('/repo'),str(path))==path


def test_alternate_config_loads(tmp_path):
    path=tmp_path/'alternate.yaml';path.write_text(yaml.safe_dump(minimal()),encoding='utf-8')
    config=load_config(path);validate_config(config)
    assert config.sequence['name']=='test'


def test_missing_config_raises(tmp_path):
    with pytest.raises(FileNotFoundError,match='config.yaml not found'):
        load_config(tmp_path/'missing.yaml')


def test_invalid_yaml_raises(tmp_path):
    path=tmp_path/'bad.yaml';path.write_text('sequence: [unterminated',encoding='utf-8')
    with pytest.raises(yaml.YAMLError): load_config(path)


@pytest.mark.parametrize('enzyme',["RNase_A","RNase_T1"])
def test_rnase_configs_load_without_changing_enzyme(tmp_path,enzyme):
    path=tmp_path/f'{enzyme}.yaml';path.write_text(yaml.safe_dump(minimal(enzyme)),encoding='utf-8')
    config=load_config(path);validate_config(config)
    assert config.digestion['enzyme']==enzyme


def test_startup_check_accepts_selected_config_path(tmp_path):
    from unittest.mock import Mock
    from rna_masshunter.startup_check import run_startup_check
    root=Path.cwd();path=tmp_path/'selected.yaml';path.write_text(yaml.safe_dump(minimal()),encoding='utf-8')
    config=load_config(path)
    config.project.update({"output_dir":str(tmp_path/'out'),"log_dir":str(tmp_path/'log'),"cache_dir":str(tmp_path/'cache')})
    warnings=[];logger=Mock()
    run_startup_check(root,config,logger,warnings,config_path=path)
    logged=[call.args for call in logger.info.call_args_list]
    assert ("Path OK: %s",path) in logged
