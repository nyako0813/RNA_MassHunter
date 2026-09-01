import builtins

from rna_masshunter.resource_utils import get_maximum_rss_mib


def test_returns_positive_float_on_current_posix_environment():
    value = get_maximum_rss_mib()
    assert isinstance(value, float)
    assert value > 0


def test_returns_none_when_resource_module_unavailable(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "resource":
            raise ImportError("simulated: resource unavailable (e.g. Windows)")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    assert get_maximum_rss_mib() is None


if __name__ == "__main__":
    test_returns_positive_float_on_current_posix_environment()
    print("resource_utils tests: OK (run under pytest for the monkeypatch case)")
