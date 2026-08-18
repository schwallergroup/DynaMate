import logging
import sys
from pathlib import Path
from src import constants
import json
import pathlib
from datetime import datetime


def get_class_logger(class_name: str, log_file: Path = None, log_to_file: bool = True) -> logging.Logger:
    """
    Create or retrieve a logger specific to a class.
    When log_file is provided all loggers in the run share the same file.
    """
    logger = logging.getLogger(class_name)
    logger.setLevel(logging.INFO)

    # Avoid adding duplicate handlers
    if not logger.handlers:
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        if log_to_file and log_file is not None:
            log_file.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

        # also print to stdout
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    return logger


def append_jsonl(data, filename):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
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
