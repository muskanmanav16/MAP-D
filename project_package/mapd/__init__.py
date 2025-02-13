"""Top-level package for mapd. can be included in init file"""


import os
import logging

from pathlib import Path
from sqlalchemy import create_engine
import pymysql
import logging
from sqlalchemy import create_engine
from pathlib import Path


HOME = Path.home()
PROJECT_DIR = HOME.joinpath(".gp2_plab2", "group2")
LOG_DIR = PROJECT_DIR.joinpath("log")
DB_PATH = PROJECT_DIR.joinpath("gp2_plab2.db")

DATA_DIR = PROJECT_DIR.joinpath("data")
PUBMED_DIR = DATA_DIR.joinpath("pubmed")
TEST_PUBMED_DIR = DATA_DIR.joinpath("pubmed_test")
TEST_DB_PATH = Path(__file__).parent.parent.joinpath("tests","data","Test_DB.db")

# Make directories
LOG_DIR.mkdir(exist_ok=True, parents=True)
PUBMED_DIR.mkdir(exist_ok=True, parents=True)
TEST_PUBMED_DIR.mkdir(exist_ok=True, parents=True)


# Logging Configuration
LOG_FILE_PATH = LOG_DIR.joinpath("group2_log.log")
logging.basicConfig(filename=LOG_FILE_PATH,
                    level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# SQLite init
CONN_STRING = f"sqlite:///{DB_PATH}"

engine = create_engine(CONN_STRING)
