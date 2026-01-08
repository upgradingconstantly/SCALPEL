"""
Base class for modality plugins.

Re-exported from core.plugin_registry for convenience.
"""

from scalpel.core.plugin_registry import ModalityPlugin, register_modality

__all__ = ["ModalityPlugin", "register_modality"]
