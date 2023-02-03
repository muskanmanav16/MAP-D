from pathlib import Path
import pytest
from mapd import DATA_DIR, LOG_DIR, LOG_FILE_PATH, DB_PATH


class TestInit:
    """Unit tests for variables in mapd-init"""

    def test_cache_dirs(self):
        """Check if the cache directories exist."""
        assert Path(DATA_DIR).is_dir()
        assert Path(LOG_DIR).is_dir()

    def test_log_file(self):
        """Check if the log file is made."""
        assert Path(LOG_FILE_PATH).is_file()

    def test_db_file(self):
        """Check if the database is created."""
        assert Path(DB_PATH).is_file()

