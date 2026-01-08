"""
Three-tier caching system for SCALPEL.

L1: In-memory LRU cache (fastest, limited size)
L2: Redis (optional, shared across processes)
L3: SQLite/DuckDB (persistent, large capacity)
"""

from __future__ import annotations

import hashlib
import json
import pickle
from functools import wraps
from typing import Any, Callable, Optional, TypeVar, ParamSpec
from collections import OrderedDict
import threading

P = ParamSpec("P")
T = TypeVar("T")


class LRUCache:
    """Thread-safe LRU cache for L1."""
    
    def __init__(self, max_size: int = 10000):
        self.max_size = max_size
        self.cache: OrderedDict[str, Any] = OrderedDict()
        self.lock = threading.Lock()
    
    def get(self, key: str) -> Optional[Any]:
        """Get value from cache, moving to end (most recent)."""
        with self.lock:
            if key in self.cache:
                self.cache.move_to_end(key)
                return self.cache[key]
            return None
    
    def set(self, key: str, value: Any) -> None:
        """Set value in cache, evicting oldest if needed."""
        with self.lock:
            if key in self.cache:
                self.cache.move_to_end(key)
            self.cache[key] = value
            while len(self.cache) > self.max_size:
                self.cache.popitem(last=False)
    
    def delete(self, key: str) -> None:
        """Remove key from cache."""
        with self.lock:
            self.cache.pop(key, None)
    
    def clear(self) -> None:
        """Clear entire cache."""
        with self.lock:
            self.cache.clear()
    
    def __len__(self) -> int:
        return len(self.cache)


class RedisCache:
    """Redis cache for L2 (optional)."""
    
    def __init__(self, url: str = "redis://localhost:6379/0", ttl: int = 3600):
        self.url = url
        self.ttl = ttl
        self._client: Optional[Any] = None
        self._available = False
        self._init_client()
    
    def _init_client(self) -> None:
        """Initialize Redis client if available."""
        try:
            import redis
            self._client = redis.from_url(self.url)
            self._client.ping()
            self._available = True
        except Exception:
            self._available = False
    
    @property
    def available(self) -> bool:
        return self._available
    
    def get(self, key: str) -> Optional[Any]:
        """Get value from Redis."""
        if not self._available or not self._client:
            return None
        try:
            data = self._client.get(f"gc:{key}")
            if data:
                return pickle.loads(data)
            return None
        except Exception:
            return None
    
    def set(self, key: str, value: Any) -> None:
        """Set value in Redis with TTL."""
        if not self._available or not self._client:
            return
        try:
            data = pickle.dumps(value)
            self._client.setex(f"gc:{key}", self.ttl, data)
        except Exception:
            pass
    
    def delete(self, key: str) -> None:
        """Remove key from Redis."""
        if not self._available or not self._client:
            return
        try:
            self._client.delete(f"gc:{key}")
        except Exception:
            pass


class DuckDBCache:
    """DuckDB cache for L3 (persistent)."""
    
    def __init__(self, db_path: str = ":memory:"):
        self.db_path = db_path
        self._connection: Optional[Any] = None
        self._init_db()
    
    def _init_db(self) -> None:
        """Initialize DuckDB and create cache table."""
        try:
            import duckdb
            self._connection = duckdb.connect(self.db_path)
            self._connection.execute("""
                CREATE TABLE IF NOT EXISTS cache (
                    key VARCHAR PRIMARY KEY,
                    value BLOB,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
        except Exception:
            self._connection = None
    
    def get(self, key: str) -> Optional[Any]:
        """Get value from DuckDB."""
        if not self._connection:
            return None
        try:
            result = self._connection.execute(
                "SELECT value FROM cache WHERE key = ?", [key]
            ).fetchone()
            if result:
                return pickle.loads(result[0])
            return None
        except Exception:
            return None
    
    def set(self, key: str, value: Any) -> None:
        """Set value in DuckDB."""
        if not self._connection:
            return
        try:
            data = pickle.dumps(value)
            self._connection.execute(
                """INSERT OR REPLACE INTO cache (key, value) VALUES (?, ?)""",
                [key, data]
            )
        except Exception:
            pass
    
    def delete(self, key: str) -> None:
        """Remove key from DuckDB."""
        if not self._connection:
            return
        try:
            self._connection.execute("DELETE FROM cache WHERE key = ?", [key])
        except Exception:
            pass


class CacheManager:
    """
    Three-tier cache manager.
    
    Automatically promotes/demotes across tiers.
    """
    
    _instance: Optional["CacheManager"] = None
    
    def __new__(cls) -> "CacheManager":
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        if self._initialized:
            return
        
        from scalpel.config import get_config
        config = get_config()
        
        # Initialize tiers
        self.l1 = LRUCache(max_size=config.cache.memory_cache_size)
        
        if config.cache.redis_enabled:
            self.l2: Optional[RedisCache] = RedisCache(
                url=config.cache.redis_url,
                ttl=config.cache.redis_ttl
            )
        else:
            self.l2 = None
        
        self.l3 = DuckDBCache(str(config.cache.db_cache_path))
        
        self._initialized = True
    
    def get(self, key: str) -> Optional[Any]:
        """Get value, checking all tiers."""
        # L1
        value = self.l1.get(key)
        if value is not None:
            return value
        
        # L2
        if self.l2 and self.l2.available:
            value = self.l2.get(key)
            if value is not None:
                self.l1.set(key, value)  # Promote to L1
                return value
        
        # L3
        value = self.l3.get(key)
        if value is not None:
            self.l1.set(key, value)  # Promote to L1
            if self.l2 and self.l2.available:
                self.l2.set(key, value)  # Promote to L2
            return value
        
        return None
    
    def set(self, key: str, value: Any, persist: bool = True) -> None:
        """Set value in all appropriate tiers."""
        self.l1.set(key, value)
        
        if self.l2 and self.l2.available:
            self.l2.set(key, value)
        
        if persist:
            self.l3.set(key, value)
    
    def delete(self, key: str) -> None:
        """Delete value from all tiers."""
        self.l1.delete(key)
        if self.l2:
            self.l2.delete(key)
        self.l3.delete(key)
    
    def clear_memory(self) -> None:
        """Clear in-memory cache only."""
        self.l1.clear()


def make_cache_key(*args: Any, **kwargs: Any) -> str:
    """Generate cache key from function arguments."""
    key_data = json.dumps({"args": args, "kwargs": kwargs}, sort_keys=True, default=str)
    return hashlib.sha256(key_data.encode()).hexdigest()[:32]


def cached(persist: bool = True) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """
    Decorator for caching function results.
    
    Args:
        persist: Whether to persist to L3 (DuckDB)
    
    Example:
        @cached()
        def expensive_computation(x, y):
            ...
    """
    def decorator(func: Callable[P, T]) -> Callable[P, T]:
        cache_manager = CacheManager()
        
        @wraps(func)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            key = f"{func.__module__}.{func.__name__}:{make_cache_key(*args, **kwargs)}"
            
            # Try cache
            result = cache_manager.get(key)
            if result is not None:
                return result
            
            # Compute and cache
            result = func(*args, **kwargs)
            cache_manager.set(key, result, persist=persist)
            return result
        
        return wrapper
    return decorator
