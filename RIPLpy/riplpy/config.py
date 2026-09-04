# -*- coding: utf-8 -*-
"""A configuration file for RIPL database access.

This file provides a centralized location for managing file paths and other
configuration settings for the RIPLpy library.

Path Resolution Priority:
    1. Programmatically set via set_path()
    2. Environment variable RIPL_LOCATION
    3. Config file ~/.riplpyrc
    4. Auto-detected common locations

Attributes:
    RIPL_PATH (str): The absolute path to the RIPL database on the local
        filesystem. Use get_path() to retrieve and set_path() to modify.
"""

# OS
import os

# Warnings
import warnings

# ========================
# Path Management
# ========================

# Internal storage for the RIPL path (set programmatically)
_ripl_path_override = None

# Common locations to check for RIPL database
_COMMON_RIPL_LOCATIONS = [
    '/usr/local/share/RIPL',
    '/usr/share/RIPL',
    '/opt/RIPL',
    os.path.expanduser('~/RIPL'),
    os.path.expanduser('~/ripl'),
    os.path.expanduser('~/Data/RIPL'),
    os.path.expanduser('~/data/RIPL'),
]

# Config file location
_CONFIG_FILE = os.path.expanduser('~/.riplpyrc')


def _read_config_file() -> str | None:
    """Read RIPL path from config file if it exists."""
    if os.path.exists(_CONFIG_FILE):
        try:
            with open(_CONFIG_FILE, 'r', encoding='utf-8', errors='replace') as f:
                for line in f:
                    line = line.strip()
                    # Skip comments and empty lines
                    if not line or line.startswith('#'):
                        continue
                    # Look for RIPL_PATH = /path/to/ripl or just /path/to/ripl
                    if '=' in line:
                        key, value = line.split('=', 1)
                        if key.strip().upper() == 'RIPL_PATH':
                            path = value.strip().strip('"').strip("'")
                            if os.path.isdir(path):
                                return path
                    elif os.path.isdir(line):
                        return line
        except (IOError, OSError):
            pass
    return None


def _auto_detect_ripl_path() -> str | None:
    """Try to auto-detect RIPL installation in common locations."""
    for location in _COMMON_RIPL_LOCATIONS:
        if os.path.isdir(location):
            # Verify it looks like a RIPL directory (has expected subdirs)
            expected_subdirs = ['masses', 'levels', 'densities']
            if all(os.path.isdir(os.path.join(location, subdir)) for subdir in expected_subdirs):
                return location
    return None


def _resolve_ripl_path() -> str:
    """Resolve RIPL path from available sources in priority order.

    Priority:
        1. Programmatically set via set_path()
        2. Environment variable RIPL_LOCATION
        3. Config file ~/.riplpyrc
        4. Auto-detected common locations
    """
    # 1. Check programmatic override
    if _ripl_path_override is not None:
        return _ripl_path_override

    # 2. Check environment variable
    env_path = os.environ.get('RIPL_LOCATION', '')
    if env_path and os.path.isdir(env_path):
        return env_path

    # 3. Check config file
    config_path = _read_config_file()
    if config_path:
        return config_path

    # 4. Try auto-detection
    auto_path = _auto_detect_ripl_path()
    if auto_path:
        return auto_path

    # Nothing found
    return ''


def set_path(path: str) -> None:
    """Set the RIPL database path programmatically.

    This takes highest priority over environment variables and config files.

    Args:
        path: Absolute path to the RIPL database directory

    Raises:
        ValueError: If the path does not exist or is not a directory

    Examples:
        >>> import riplpy
        >>> riplpy.config.set_path('/path/to/RIPL-4')
        >>> riplpy.load()
    """
    global _ripl_path_override, RIPL_PATH

    if not os.path.isdir(path):
        raise ValueError(f"Path does not exist or is not a directory: {path}")

    _ripl_path_override = os.path.abspath(path)
    RIPL_PATH = _ripl_path_override


def get_path() -> str:
    """Get the current RIPL database path.

    Returns:
        The resolved RIPL path, or empty string if not configured

    Examples:
        >>> import riplpy
        >>> print(riplpy.config.get_path())
        '/usr/local/share/RIPL-4'
    """
    return _resolve_ripl_path()


def resolve_directory(directory: str = None) -> str:
    """Resolve a section loader's ``directory`` argument.

    If ``directory`` is provided it is returned unchanged. Otherwise the
    configured RIPL path is used (set_path(), RIPL_LOCATION env var,
    ~/.riplpyrc, or an auto-detected location). This gives every section
    ``load()`` the same fallback behaviour as ``riplpy.load()``.

    Args:
        directory: Explicit path, or None to use the configured path.

    Returns:
        The resolved RIPL database directory.

    Raises:
        ValueError: If no directory is given and no path is configured.
    """
    if directory is not None:
        return directory
    resolved = _resolve_ripl_path()
    if not resolved:
        raise ValueError(
            "No RIPL directory specified and no path configured. "
            "Either pass a directory to load(), or configure the path first:\n"
            "  - riplpy.set_path('/path/to/RIPL')\n"
            "  - Set RIPL_LOCATION environment variable\n"
            "  - Create ~/.riplpyrc config file"
        )
    return resolved


def clear_path() -> None:
    """Clear the programmatically set path, reverting to other sources."""
    global _ripl_path_override, RIPL_PATH
    _ripl_path_override = None
    RIPL_PATH = _resolve_ripl_path()


def write_config_file(path: str) -> None:
    """Write the RIPL path to the config file for persistence.

    Args:
        path: Path to write to the config file

    Examples:
        >>> import riplpy
        >>> riplpy.config.write_config_file('/path/to/RIPL-4')
        # Creates ~/.riplpyrc with the path
    """
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        raise ValueError(f"Path does not exist or is not a directory: {path}")

    with open(_CONFIG_FILE, 'w', encoding='utf-8') as f:
        f.write(f"# RIPLpy configuration file\n")
        f.write(f"# Generated automatically\n")
        f.write(f"RIPL_PATH = {path}\n")


# Initialize RIPL_PATH
# ========================
# Resolve the path on module load
RIPL_PATH = _resolve_ripl_path()

# Only warn if no path could be resolved
if not RIPL_PATH:
    warnings.warn(
        "RIPL database path not configured. Options to set it:\n"
        "  1. riplpy.set_path('/path/to/RIPL')\n"
        "  2. Set RIPL_LOCATION environment variable\n"
        "  3. Create ~/.riplpyrc with: RIPL_PATH = /path/to/RIPL\n"
        "  4. Install RIPL to a standard location (~/RIPL, /usr/local/share/RIPL)",
        UserWarning
    )

# Standardized Relative Paths for ASCII Files

# This section provides a centralized mapping of descriptive names to the
# relative paths of the ASCII data files within the RIPL database. This
# helps to avoid hardcoding paths throughout the codebase.

DATA_FILES = {

    None: "", # Default key for no file path

    # Densities section
    # =================

    # Level densities
    'densities': os.path.join('densities', 'data', 'densities.dat'),
    # Note: BSFG (level-densities-bfmeff.dat) and CT (level-densities-ctmeff.dat)
    # are not present in the public RIPL-4 github release. The keys are retained
    # so users with the older full release can still load these databases.
    'densities_bsfg': os.path.join('densities', 'level-densities-bfmeff.dat'),
    'densities_ct': os.path.join('densities', 'level-densities-ctmeff.dat'),
    'densities_egsm': os.path.join('densities', 'total', 'level-densities-egsm.dat'),
    'densities_egsm_norm': os.path.join('densities', 'total', 'level-densities-egsm-norm.dat'),

    # Microscopic combinatorial level density tables (per-Z .tab files)
    # under densities/total/<model>/zXXX.tab
    'densities_bsk14_comb_dir': os.path.join('densities', 'total', 'bsk14-comb'),
    'densities_bskg3_comb_dir': os.path.join('densities', 'total', 'bskg3-comb'),
    'densities_qrpabe_dir': os.path.join('densities', 'total', 'qrpabe'),
    'densities_thfb_comb_dir': os.path.join('densities', 'total', 'thfb-comb'),

    # Fission section
    # ===============

    # Fission Barriers
    'fission_barriers_empirical': os.path.join('fission', 'empirical-barriers-ripl4.dat'),
    'fission_barriers_empirical_new': os.path.join('fission', 'empirical-barriers-ripl4.dat'),
    'fission_barriers_empire': os.path.join('fission', 'empirical-barriers-new-EMPIRE.dat'),
    'fission_barriers_hfb': os.path.join('fission', 'empirical-hfb-barriers.dat'),
    'fission_barriers_bskg3': os.path.join('fission', 'barriers-bskg3.dat'),
    'fission_barriers_d1m': os.path.join('fission', 'barriers-d1m_lep.dat'),

    # Fission paths (per-isotope directories)
    'fission_path_bskg3_dir': os.path.join('fission', 'hfbpath-bskg3'),
    'fission_path_d1m_dir': os.path.join('fission', 'hfbpath-d1m'),

    # RMF fission paths (per-nucleus files in axial/triaxial directories)
    'fission_rmf_axial_dir': os.path.join('fission', 'RMF', 'Path_Axial'),
    'fission_rmf_triaxial_dir': os.path.join('fission', 'RMF', 'Path_Triaxial'),

    # Fission saddle-point level densities (per-Z files in subdirs Max1/Max2/Max3/Min1/Min2)
    'fission_nld_bskg3_dir': os.path.join('fission', 'nld-fis-bskg3'),

    # Gamma section
    # =============
    #
    # The RIPL-4 github release reorganises the gamma layout substantially.
    # Legacy single-file products (e.g. ``gamma/gamma-strength-micro``,
    # ``gamma/gdr-parameters-theor.dat``, and ``gamma/data/...``) are no longer
    # shipped. The keys below point at the new locations; absent legacy keys are
    # retained for backwards compatibility but the corresponding loaders are
    # tolerant to the underlying file being missing.

    # Gamma-ray strength functions (legacy single-file product; not in github)
    'gsf_data': os.path.join('gamma', 'data', 'gsf', 'gsf_Z{Z:03d}A{A:03d}.dat'),

    # Theoretical / experimental GDR parameter compilations
    'gdr_parameters_theor': os.path.join('gamma', 'gdr-parameters-theor.dat'),
    'gdr_parameters_exp': os.path.join('gamma', 'data', 'ripl3-gdr-parameters&errors-exp-SLO.dat'),

    # RIPL-4 recommended experimental SLO / SMLO GDR parameters
    'gdr_parameters_recommended_exp_slo': os.path.join(
        'gamma', 'gdr_parameters_exp_new', 'gdr_parameters_recommended_exp_slo.dat'
    ),
    'gdr_parameters_recommended_exp_smlo': os.path.join(
        'gamma', 'gdr_parameters_exp_new', 'gdr_parameters_recommended_exp_smlo.dat'
    ),
    # Experimental values WITH errors (two-line records)
    'gdr_parameters_errors_exp_slo': os.path.join(
        'gamma', 'gdr_parameters_exp_new', 'gdr_parameters&errors_exp_slo.dat'
    ),
    'gdr_parameters_errors_exp_smlo': os.path.join(
        'gamma', 'gdr_parameters_exp_new', 'gdr_parameters&errors_exp_smlo.dat'
    ),
    # Broader SLO / SMLO fits including systematics (~8980 entries)
    'gdr_parameters_systematics_slo': os.path.join(
        'gamma', 'gdr_parameters_exp&systematics', 'gdr-parameters_exp&systematics_slo.dat'
    ),
    'gdr_parameters_systematics_smlo': os.path.join(
        'gamma', 'gdr_parameters_exp&systematics', 'gdr-parameters_exp&systematics_smlo.dat'
    ),

    # D1M+QRPA per-nucleus gamma strength predictions
    'gsf_d1m_dir': os.path.join('gamma', 'd1m'),

    # SMLO E1 photoabsorption tables (per-nucleus, fe1_the_<Z>_<A>_photoabs_h_SMLO.dat)
    'gsf_smlo_e1_dir': os.path.join('gamma', 'smlo_E1'),
    # SMLO M1 strength tables (per-Z files, e.g. z026_m1)
    'gsf_smlo_m1_dir': os.path.join('gamma', 'smlo_M1'),
    # TLO data (per-Z FE1_zXXX.dat files)
    'gsf_tlo_dir': os.path.join('gamma', 'tlo'),

    # Photon Strength Function Database (PSFDatabase-v2024.1, multi-format)
    'psf_database_dir': os.path.join('gamma', 'PSFDatabase-v2024.1'),

    # Levels section
    # ==============

    # Levels
    'levels_param': os.path.join('levels', 'levels-param.data'),
    'discrete_levels': os.path.join('levels', 'z{Z:03d}.dat'),

    # Masses section
    # ==============

    # Mass models
    'mass_ame20': os.path.join('masses', 'mass-ame20.dat'),
    'mass_bskg3': os.path.join('masses', 'mass-bskg3.dat'),
    'mass_d1m': os.path.join('masses', 'mass-d1m.dat'),
    'mass_frdm12': os.path.join('masses', 'mass-frdm12.dat'),
    'mass_frdm95': os.path.join('masses', 'mass-frdm95.dat'),
    'mass_hfb14': os.path.join('masses', 'mass-hfb14.dat'),
    'mass_hfb27': os.path.join('masses', 'mass-hfb27.dat'),
    'mass_ws4': os.path.join('masses', 'mass-ws4.dat'),

    # Deformations
    'deformations_gs_exp': os.path.join('masses', 'gs-deformations-exp.dat'),

    # Abundances
    'abundances': os.path.join('masses', 'abundance.dat'),

    # Optical Model section
    # =====================

    # Optical Model Parameters
    'om_deformations': os.path.join('optical', 'data', 'om-deformations.dat'),
    'om_parameter_u': os.path.join('optical', 'data', 'om-parameter-u.dat'),

    # Resonances section
    # ==================

    # Resonances
    'resonances_swave': os.path.join('resonances', 'resonances_L0.dat'),
    'resonances_pwave': os.path.join('resonances', 'resonances_L1.dat'),
}

def get_data_file_path(key: str, **kwargs) -> str:
    """Constructs the full path to a data file.

    Arguments:
        key (str): The key corresponding to the data file in the DATA_FILES dict.
        **kwargs: Formatting arguments for paths that require them (e.g., Z and A).

    Returns:
        str: The absolute path to the data file.
    """
    if key not in DATA_FILES:
        raise KeyError(f"Data file key '{key}' not found in configuration.")
    
    relative_path = DATA_FILES[key]
    
    # Format the path if needed
    if '{' in relative_path:
        relative_path = relative_path.format(**kwargs)
        
    if not RIPL_PATH:
        # Return relative path if RIPL_PATH is not set, which might still work
        # if the application is run from a suitable directory.
        return relative_path

    return os.path.join(RIPL_PATH, relative_path)

