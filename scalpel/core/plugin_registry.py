"""
Plugin registry for modality designers.

Allows dynamic registration and discovery of editing modality plugins.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, Type, List, Optional, Any, Callable
from functools import wraps

from scalpel.models.data_classes import (
    ResolvedTarget,
    DesignedGuide,
    DesignResults,
    EfficiencyScore,
)
from scalpel.models.enums import EditModality, CasVariantType


class ModalityPlugin(ABC):
    """
    Base class for modality-specific designers.
    
    Each modality (knockout, CRISPRi, base editing, etc.) implements
    this interface to provide modality-specific guide design logic.
    """
    
    # Class attributes to be set by subclasses
    modality: EditModality
    display_name: str
    description: str
    supported_cas_variants: List[CasVariantType] = [CasVariantType.SPCAS9]
    
    # Scoring weight profile (modality-specific)
    scoring_weights: Dict[str, float] = {}
    
    @abstractmethod
    def design(
        self,
        target: ResolvedTarget,
        n_guides: int = 10,
        **kwargs: Any,
    ) -> DesignResults:
        """
        Design guides for the given target.
        
        Args:
            target: Resolved target specification
            n_guides: Number of top guides to return
            **kwargs: Modality-specific parameters
        
        Returns:
            DesignResults with ranked guides
        """
        pass
    
    @abstractmethod
    def score_position(
        self,
        guide: DesignedGuide,
        target: ResolvedTarget,
    ) -> float:
        """
        Score guide position for this modality.
        
        Different modalities have different position preferences:
        - Knockout: Early in CDS, avoid last exon
        - CRISPRi: Near TSS
        - CRISPRa: Upstream of TSS
        - Base editing: Edit in optimal window
        
        Returns:
            Position score in [0, 1]
        """
        pass
    
    def validate_target(self, target: ResolvedTarget) -> List[str]:
        """
        Validate target is suitable for this modality.
        
        Returns:
            List of validation errors (empty if valid)
        """
        errors = []
        
        # Check cas variant compatibility
        cas_variant = target.specification.cas_variant
        if cas_variant not in self.supported_cas_variants:
            errors.append(
                f"{self.display_name} does not support {cas_variant.value}. "
                f"Supported: {[cv.value for cv in self.supported_cas_variants]}"
            )
        
        return errors
    
    def get_default_parameters(self) -> Dict[str, Any]:
        """Get default parameters for this modality."""
        return {}


class ModalityRegistry:
    """
    Registry for modality plugins.
    
    Enables dynamic discovery and instantiation of modality designers.
    """
    
    _instance: Optional["ModalityRegistry"] = None
    _plugins: Dict[EditModality, Type[ModalityPlugin]] = {}
    
    def __new__(cls) -> "ModalityRegistry":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    @classmethod
    def register(cls, modality: EditModality) -> Callable[[Type[ModalityPlugin]], Type[ModalityPlugin]]:
        """
        Decorator to register a modality plugin.
        
        Example:
            @ModalityRegistry.register(EditModality.KNOCKOUT)
            class KnockoutDesigner(ModalityPlugin):
                ...
        """
        def decorator(plugin_class: Type[ModalityPlugin]) -> Type[ModalityPlugin]:
            plugin_class.modality = modality
            cls._plugins[modality] = plugin_class
            return plugin_class
        return decorator
    
    @classmethod
    def get(cls, modality: EditModality) -> Optional[Type[ModalityPlugin]]:
        """Get plugin class for a modality."""
        return cls._plugins.get(modality)
    
    @classmethod
    def get_instance(cls, modality: EditModality) -> Optional[ModalityPlugin]:
        """Get instantiated plugin for a modality."""
        plugin_class = cls.get(modality)
        if plugin_class:
            return plugin_class()
        return None
    
    @classmethod
    def list_modalities(cls) -> List[EditModality]:
        """List all registered modalities."""
        return list(cls._plugins.keys())
    
    @classmethod
    def list_plugins(cls) -> Dict[EditModality, Type[ModalityPlugin]]:
        """Get all registered plugins."""
        return cls._plugins.copy()


# Convenience decorator
def register_modality(modality: EditModality) -> Callable[[Type[ModalityPlugin]], Type[ModalityPlugin]]:
    """
    Convenience decorator for registering modality plugins.
    
    Example:
        @register_modality(EditModality.KNOCKOUT)
        class KnockoutDesigner(ModalityPlugin):
            display_name = "CRISPR Knockout"
            ...
    """
    return ModalityRegistry.register(modality)
