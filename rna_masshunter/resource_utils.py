"""Cross-platform process resource usage utilities.

The stdlib ``resource`` module is POSIX-only (Linux/macOS) and does not
exist on Windows. This module centralizes the one RSS-memory lookup used
across the codebase so callers never import ``resource`` directly.
"""
from __future__ import annotations


def get_maximum_rss_mib() -> float | None:
    """Return the process's peak resident set size (RSS) in MiB.

    Returns ``None`` on platforms where the POSIX ``resource`` module is
    unavailable (currently: Windows). This is a diagnostic-only metric;
    callers must handle ``None`` explicitly rather than assuming a float.
    """
    try:
        import resource
    except ImportError:
        return None
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
