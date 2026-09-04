
# ========================

# Functools
from functools import partial as _partial

# RIPLpy
import riplpy.config as _config
import riplpy.collections as _c
import riplpy.db as _db
from . import pwave as _pwave

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'resonances_swave'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

# The s- and p-wave files share the same RIPL-4 format (Goriely 2025), so
# the reader, writer, and Entry dataclass are imported directly from pwave.
FORTRAN_FORMAT = _pwave.FORTRAN_FORMAT
FILE_HEADER    = _pwave.FILE_HEADER

read_ascii_file  = _pwave.read_ascii_file
write_ascii_file = _pwave.write_ascii_file

Entry = _pwave.Entry


class Database(_db.NuclideDatabase):

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
