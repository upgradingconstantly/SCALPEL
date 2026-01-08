"""Modality plugins package."""

# Import all modality plugins to trigger registration
from scalpel.modalities.base import ModalityPlugin
from scalpel.modalities.knockout import KnockoutDesigner
from scalpel.modalities.interference import InterferenceDesigner
from scalpel.modalities.activation import ActivationDesigner
from scalpel.modalities.base_editor import BaseEditorDesigner
from scalpel.modalities.prime_editor import PrimeEditorDesigner

__all__ = [
    "ModalityPlugin",
    "KnockoutDesigner",
    "InterferenceDesigner",
    "ActivationDesigner",
    "BaseEditorDesigner",
    "PrimeEditorDesigner",
]
