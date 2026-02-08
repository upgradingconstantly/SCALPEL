"""
CRISPRon-style efficiency prediction model.

Uses a CNN architecture similar to CRISPRon for gRNA efficiency prediction.
Can load pre-trained weights or use rule-based fallback.

References:
- CRISPRon: https://github.com/RTH-tools/crispron
- Xiang et al. 2021 - Data integration and deep learning
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, List, Tuple
import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

logger = logging.getLogger(__name__)


# One-hot encoding for DNA sequences
NUCLEOTIDE_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}


def one_hot_encode(sequence: str) -> np.ndarray:
    """
    One-hot encode a DNA sequence.
    
    Args:
        sequence: DNA sequence (A, C, G, T, N)
    
    Returns:
        numpy array of shape (4, len(sequence))
    """
    seq = sequence.upper()
    encoding = np.zeros((4, len(seq)), dtype=np.float32)
    for i, base in enumerate(seq):
        idx = NUCLEOTIDE_MAP.get(base, 0)
        encoding[idx, i] = 1.0
    return encoding


class CRISPRonCNN(nn.Module):
    """
    CNN architecture for gRNA efficiency prediction.
    
    Based on CRISPRon architecture with:
    - Convolutional layers to capture sequence motifs
    - Global pooling for position-invariant features
    - Dense layers for final prediction
    """
    
    def __init__(self, sequence_length: int = 30):
        super().__init__()
        
        # Convolutional layers
        self.conv1 = nn.Conv1d(4, 64, kernel_size=5, padding=2)
        self.conv2 = nn.Conv1d(64, 128, kernel_size=5, padding=2)
        self.conv3 = nn.Conv1d(128, 64, kernel_size=3, padding=1)
        
        # Batch normalization
        self.bn1 = nn.BatchNorm1d(64)
        self.bn2 = nn.BatchNorm1d(128)
        self.bn3 = nn.BatchNorm1d(64)
        
        # Dropout
        self.dropout = nn.Dropout(0.3)
        
        # Dense layers
        self.fc1 = nn.Linear(64, 64)
        self.fc2 = nn.Linear(64, 1)
    
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Forward pass.
        
        Args:
            x: Input tensor of shape (batch, 4, sequence_length)
        
        Returns:
            Prediction tensor of shape (batch, 1)
        """
        # Convolutional blocks
        x = F.relu(self.bn1(self.conv1(x)))
        x = F.relu(self.bn2(self.conv2(x)))
        x = F.relu(self.bn3(self.conv3(x)))
        
        # Global average pooling
        x = x.mean(dim=2)  # (batch, 64)
        
        # Dense layers
        x = self.dropout(x)
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = torch.sigmoid(self.fc2(x))
        
        return x


class CRISPRonPredictor:
    """
    CRISPRon-style efficiency predictor.
    
    Loads pre-trained weights if available, otherwise uses rule-based scoring
    with optional training capability.
    """
    
    # Default model path
    DEFAULT_MODEL_PATH = Path.home() / ".scalpel" / "models" / "crispron" / "model.pt"
    
    def __init__(
        self,
        model_path: Optional[Path] = None,
        use_rule_fallback: bool = True,
    ):
        """
        Initialize the predictor.
        
        Args:
            model_path: Path to pre-trained model weights
            use_rule_fallback: If True, use rule-based scoring when model unavailable
        """
        self.model_path = model_path or self.DEFAULT_MODEL_PATH
        self.use_rule_fallback = use_rule_fallback
        self.model: Optional[CRISPRonCNN] = None
        self.device = "cpu"
        
        # Try to load model
        self._load_model()
    
    def _load_model(self) -> None:
        """Load the pre-trained model if available."""
        if not TORCH_AVAILABLE:
            logger.warning("PyTorch not available, using rule-based fallback")
            return
        
        if self.model_path.exists():
            # B4.2: Validate checksum if manifest exists
            if not self._validate_checksum():
                logger.warning("Model checksum validation failed, using rule-based fallback")
                return
            
            try:
                self.model = CRISPRonCNN()
                state_dict = torch.load(self.model_path, map_location="cpu")
                self.model.load_state_dict(state_dict)
                self.model.eval()
                
                # Use GPU if available
                if torch.cuda.is_available():
                    self.device = "cuda"
                    self.model = self.model.to(self.device)
                
                logger.info(f"Loaded CRISPRon model from {self.model_path}")
            except Exception as e:
                logger.warning(f"Failed to load model: {e}")
                self.model = None
        else:
            logger.info(
                f"No pre-trained model found at {self.model_path}. "
                "Using rule-based scoring. To enable ML predictions, "
                "download weights or train the model."
            )
    
    def _validate_checksum(self) -> bool:
        """
        B4.2: Validate model checksum against manifest.
        
        Returns True if no manifest exists (assume valid) or checksum matches.
        Returns False if manifest exists but checksum doesn't match.
        """
        import hashlib
        
        manifest_path = self.model_path.parent / "manifest.json"
        if not manifest_path.exists():
            # No manifest = no validation (backwards compatible)
            return True
        
        try:
            import json
            with open(manifest_path) as f:
                manifest = json.load(f)
            
            expected_checksum = manifest.get("sha256")
            expected_version = manifest.get("version")
            
            if expected_checksum:
                # Calculate actual checksum
                sha256 = hashlib.sha256()
                with open(self.model_path, "rb") as f:
                    for chunk in iter(lambda: f.read(8192), b""):
                        sha256.update(chunk)
                actual_checksum = sha256.hexdigest()
                
                if actual_checksum != expected_checksum:
                    logger.error(
                        f"Model checksum mismatch! "
                        f"Expected: {expected_checksum[:16]}..., "
                        f"Got: {actual_checksum[:16]}..."
                    )
                    return False
                
                logger.info(f"Model checksum validated (version: {expected_version or 'unknown'})")
            
            return True
            
        except Exception as e:
            logger.warning(f"Could not validate checksum: {e}")
            return True  # Allow loading if validation code fails
    
    @property
    def is_model_available(self) -> bool:
        """Check if ML model is loaded and available."""
        return self.model is not None
    
    def predict(self, sequence: str) -> float:
        """
        Predict efficiency score for a single sequence.
        
        Args:
            sequence: 30bp sequence (20bp spacer + 3bp PAM + 7bp context)
                     or 23bp sequence (20bp spacer + 3bp PAM)
        
        Returns:
            Efficiency score between 0 and 1
        """
        # Pad to 30bp if needed
        if len(sequence) < 30:
            sequence = sequence + "N" * (30 - len(sequence))
        elif len(sequence) > 30:
            sequence = sequence[:30]
        
        if self.model is not None:
            return self._predict_ml(sequence)
        elif self.use_rule_fallback:
            return self._predict_rules(sequence)
        else:
            return 0.5  # Neutral score
    
    def predict_batch(self, sequences: List[str]) -> List[float]:
        """
        Predict efficiency scores for multiple sequences.
        
        More efficient than calling predict() repeatedly.
        
        Args:
            sequences: List of sequences
        
        Returns:
            List of efficiency scores
        """
        if not sequences:
            return []
        
        if self.model is not None:
            return self._predict_ml_batch(sequences)
        else:
            return [self.predict(seq) for seq in sequences]
    
    def _predict_ml(self, sequence: str) -> float:
        """ML prediction for single sequence."""
        return self._predict_ml_batch([sequence])[0]
    
    def _predict_ml_batch(self, sequences: List[str]) -> List[float]:
        """ML prediction for batch of sequences."""
        if self.model is None:
            raise RuntimeError("Model not loaded")
        
        # Prepare input
        padded = []
        for seq in sequences:
            if len(seq) < 30:
                seq = seq + "N" * (30 - len(seq))
            elif len(seq) > 30:
                seq = seq[:30]
            padded.append(seq)
        
        # One-hot encode
        encodings = np.stack([one_hot_encode(seq) for seq in padded])
        x = torch.tensor(encodings, dtype=torch.float32).to(self.device)
        
        # Predict
        with torch.no_grad():
            predictions = self.model(x).cpu().numpy().flatten()
        
        return predictions.tolist()
    
    def _predict_rules(self, sequence: str) -> float:
        """
        Rule-based efficiency prediction (Doench-inspired rules).
        
        Uses validated sequence features correlated with efficiency:
        - GC content in optimal range (40-70%)
        - Position-specific nucleotide preferences
        - Avoidance of Pol III terminators
        - Seed region GC content
        """
        score = 0.5  # Start neutral
        seq = sequence.upper()
        spacer = seq[:20] if len(seq) >= 20 else seq
        
        # GC content scoring
        gc_count = spacer.count("G") + spacer.count("C")
        gc_fraction = gc_count / len(spacer)
        
        if 0.40 <= gc_fraction <= 0.70:
            score += 0.15  # Optimal range
        elif gc_fraction < 0.30 or gc_fraction > 0.80:
            score -= 0.15  # Poor range
        else:
            score += 0.05  # Acceptable
        
        # Seed GC content (positions 1-12)
        seed = spacer[:12]
        seed_gc = (seed.count("G") + seed.count("C")) / len(seed)
        if 0.33 <= seed_gc <= 0.67:
            score += 0.05
        
        # Position-specific preferences (Doench Rule Set 2)
        # Position 20 (PAM-proximal): G or C preferred
        if len(spacer) >= 20:
            if spacer[19] in "GC":
                score += 0.05
        
        # Position 1: G preferred
        if spacer[0] == "G":
            score += 0.03
        
        # Penalize Pol III terminator (TTTT)
        if "TTTT" in spacer:
            score -= 0.20
        
        # Penalize G-quadruplex (GGGG)
        if "GGGG" in spacer:
            score -= 0.10
        
        # Clamp to valid range
        return max(0.0, min(1.0, score))
    
    def save_model(self, path: Optional[Path] = None) -> None:
        """Save the current model weights."""
        if self.model is None:
            raise RuntimeError("No model to save")
        
        save_path = path or self.model_path
        save_path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(self.model.state_dict(), save_path)
        logger.info(f"Saved model to {save_path}")


# Global predictor instance (lazy loaded)
_predictor: Optional[CRISPRonPredictor] = None


def get_predictor() -> CRISPRonPredictor:
    """Get the global CRISPRon predictor instance."""
    global _predictor
    if _predictor is None:
        _predictor = CRISPRonPredictor()
    return _predictor


def predict_efficiency(sequence: str) -> float:
    """
    Convenience function to predict gRNA efficiency.
    
    Args:
        sequence: 20-30bp sequence (spacer + optional PAM + context)
    
    Returns:
        Efficiency score (0-1)
    """
    return get_predictor().predict(sequence)


def predict_efficiency_batch(sequences: List[str]) -> List[float]:
    """
    Convenience function to predict efficiency for multiple sequences.
    
    Args:
        sequences: List of sequences
    
    Returns:
        List of efficiency scores
    """
    return get_predictor().predict_batch(sequences)
