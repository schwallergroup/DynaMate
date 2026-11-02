from pathlib import Path

MAX_ITERATIONS = 25
MAX_CHARACTERS_TO_LOG = 5000
SUMMARY_OUTPUT_TOKENS = 6000
MAX_CONTEXT_TOKENS = 32000
PAPER_DIR = Path(__file__).resolve().parent.parent / "my_papers"
MODEL_NAME = "openrouter/openai/gpt-5.5"
TEMPERATURE = 0.1

# Ionic environment chosen by the PrepAgent -> NaCl concentration (M) given to gmx genion (-neutral is always applied)
IONIC_ENVIRONMENTS = {
    "physiological": 0.15,  # physiological salt, e.g. human/mammalian proteins in cells or blood
    "neutral": 0.0,         # counter-ions only, no added salt, e.g. environmental proteins
}
DEFAULT_IONIC_ENVIRONMENT = "physiological"

# Equilibration (NVT then NPT) and its checks
EQUIL_TIME_NVT_PS = 100                         # NVT and NPT length
EQUIL_TIME_NPT_PS = 500                         # NVT and NPT length
EQUIL_TEST_TIME_PS = 10                     # NVT and NPT length with --test (no extension)
EQUIL_EXTENSION_PS = 100                    # added to a stage that fails its check
MAX_NVT_EXTENSIONS = 1                      # 100 ps -> at most 200 ps
MAX_NPT_EXTENSIONS = 1
NVT_TEMP_TOLERANCE_K = 5.0                  # |mean T - ref T| over the second half of NVT
REF_PRESSURE_BAR = 1.0                      # passed if |mean P - ref P| <= std P over the second half of NPT
DENSITY_RANGE_KG_M3 = (900.0, 1100.0)       # mean density over the second half of NPT
DENSITY_MAX_DRIFT_KG_M3 = 10.0              # linear drift of the density over the second half of NPT

SCRIPTS_DIR = Path(__file__).resolve().parent / "scripts"
MDP_FILES = SCRIPTS_DIR / "mdp_files"

ENV_FILE = Path(__file__).resolve().parent.parent / ".env"
DATA_DIR = Path(__file__).resolve().parent.parent / "sandbox"

AGENT_LOGS = Path(__file__).resolve().parent.parent / "agent_logs"
JSON_LOG_FILE = AGENT_LOGS / "agent_runs.jsonl"

MMPBSA_ENV_DIR = Path("/path/to/miniconda3/envs/gmxMMPBSA/bin/gmx_MMPBSA")
