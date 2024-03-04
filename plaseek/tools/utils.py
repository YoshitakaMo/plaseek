import logging
from logging import StreamHandler, Formatter
from pathlib import Path


def setup_logging(log_file: Path, mode: str = "w") -> None:
    log_file.parent.mkdir(exist_ok=True, parents=True)
    logger = logging.getLogger(__name__)
    if logger.handlers:
        for handler in logger.handlers:
            handler.close()
            logger.removeHandler(handler)
    logger.setLevel(logging.INFO)
    stream_handler = StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    handler_format = Formatter(
        "%(asctime)s %(filename)s:%(lineno)d - %(levelname)s - %(message)s"
    )
    stream_handler.setFormatter(handler_format)
    logger.addHandler(stream_handler)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(filename)s:%(lineno)d - %(levelname)s - %(message)s",
        handlers=[stream_handler, logging.FileHandler(log_file, mode=mode)],
    )
