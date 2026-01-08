"""Core utilities package."""

from scalpel.core.cache import CacheManager, cached
from scalpel.core.plugin_registry import (
    ModalityPlugin,
    ModalityRegistry,
    register_modality,
)
from scalpel.core.red_flags import (
    RedFlagDetector,
    RedFlagType,
    detect_red_flags,
    summarize_red_flags,
)

__all__ = [
    "CacheManager",
    "cached",
    "ModalityPlugin",
    "ModalityRegistry",
    "register_modality",
    "RedFlagDetector",
    "RedFlagType",
    "detect_red_flags",
    "summarize_red_flags",
]
