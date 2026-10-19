# SearchFileName.py
import glob
import os

def Search_File(pattern: str) -> str:
    """
    Search for files matching 'pattern' and return the first match.
    Example: Search_File("data/*.txt")
    """
    pattern = pattern + '*'
    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(f"No files found for pattern: {pattern}")
    # Return absolute path of the first match
    return os.path.abspath(matches[0])
