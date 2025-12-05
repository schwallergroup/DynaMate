import logging
import sys
from pathlib import Path
from src import constants
import json
import pathlib
from datetime import datetime


def get_class_logger(class_name: str, log_dir: Path = None) -> logging.Logger:
    """
    Create or retrieve a logger specific to a class.
    Each class writes to its own log file inside agent_logs/.
    """
    if log_dir is None:
        log_dir = Path(__file__).resolve().parent.parent / "agent_logs"
    log_dir.mkdir(exist_ok=True)

    log_file = log_dir / f"{class_name}.log"

    logger = logging.getLogger(class_name)
    logger.setLevel(logging.INFO)

    # Avoid adding duplicate handlers
    if not logger.handlers:
        file_handler = logging.FileHandler(log_file)
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        # also print to stdout
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    return logger


def append_jsonl(data, filename):
    """Append one JSON record per line."""
    with open(filename, "a", encoding="utf-8") as f:
        f.write(json.dumps(data, ensure_ascii=False) + "\n")


def truncate_string(string):
    if not string:
        return ""
    
    if len(string) <= 2 * constants.MAX_CHARACTERS_TO_LOG:
        return string
    
    return f"{string[:constants.MAX_CHARACTERS_TO_LOG]}... truncated ... {string[-constants.MAX_CHARACTERS_TO_LOG:]}"


def is_path_child_dir(potential_child_dir: str | Path, dir: str | Path) -> bool:
    """Ensure that the requested path stays within the sandbox."""
    if isinstance(potential_child_dir, str):
        potential_child_dir = pathlib.Path(potential_child_dir)

    if isinstance(dir, str):
        dir = pathlib.Path(dir)

    abs_potential_child = potential_child_dir.resolve()
    abs_dir = dir.resolve()

    return abs_potential_child.is_relative_to(abs_dir)


def time_now(time_format: str = "%Y%m%d_%H%M%S"):
    return datetime.now().strftime(time_format)
