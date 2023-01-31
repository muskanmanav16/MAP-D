import scispacy
import os
import sqlalchemy
import logging
from sqlalchemy import create_engine, Column, Integer, String, Table, MetaData,ForeignKey, DATE,text,select, inspect,func
from typing import Optional, Union
import en_ner_bionlp13cg_md    #The model we are going to use
import spacy
from sqlalchemy.orm import declarative_base, sessionmaker
from mapd.models import Base,Abstract,Entity
import tensorflow as tf
from sqlalchemy_utils import database_exists

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#check the device name and use GPU if it is available
device_name = tf.test.gpu_device_name()
CONN_STRING = "sqlite:///abstract.db"
from sqlalchemy.orm import sessionmaker
if device_name != '/device:GPU:0':
    print('GPU device not found. Using CPU instead.')
    tf.config.set_visible_devices([], 'GPU')
else:
    print('Found GPU at: {}'.format(device_name))

#DB connection
engine = create_engine('sqlite:///abstract.db')
Session = sessionmaker(bind=engine)
session = Session()
Base.metadata.create_all(engine)

# Load the sciSpacy model
model_name="en_ner_bionlp13cg_md"
nlp = spacy.load(model_name)

class Database:
    """For interfacing with the relational database."""

    def __init__(self, db_engine=engine):
        self.engine = db_engine
        self.session = Session(bind=self.engine)

        # When DB doesn't exist
        if not database_exists(self.engine.url):
            self.build_database()

        # When DB exists, but tables not made yet
        tables = inspect(self.engine).get_table_names()
        if not tables:
            self.build_database()

    def rebuild_database(self) -> None:
        """Burn everything and builds the database."""
        self.drop_database()
        self.build_database()

    def build_database(self) -> None:
        """Build the tables of the database."""
        logger.info("Building database...")
        Base.metadata.create_all(bind=self.engine)

    def drop_database(self) -> None:
        """Drop all of the associated tables in the database."""
        logger.warning("Dropping database...")
        Base.metadata.drop_all(bind=self.engine)
Database()
# Define a class for entity prediction

class EntityPrediction:
    def __init__(self, session):
        self.session = session

    def predict_entities(self, abstract_text):
        doc = nlp(abstract_text)
        entities = []
        for ent in doc.ents:
            entities.append({"entity": ent.text, "labels": ent.label_})
        return entities

    def insert_entities(self, abstract_id, entities):
        for entity in entities:
            e = Entity(entity=entity["entity"], labels=entity["labels"], pubmed_id=abstract_id)
            self.session.add(e)
        self.session.commit()

def get_entity_data(session):
    entity_predictor = EntityPrediction(session)
    session.execute(text('DELETE FROM Abstract WHERE id NOT IN (SELECT min(id) FROM Abstract GROUP BY pubmed_id)'))
    session.commit()
    # Query the abstract table and predict entities
    abstracts = session.query(Abstract).all()
    for abstract in abstracts:
        entities = entity_predictor.predict_entities(abstract.abstract_text)
        entity_predictor.insert_entities(abstract.id, entities)

    session.close()

def delete_duplicates_entity_table(engine):
    # Create a session
    Session = sessionmaker(bind=engine)
    session = Session()

    # Query for duplicates
    query = session.query(Entity).group_by(Entity.pubmed_id).having(func.count(Entity.id) > 1)

    # Delete the duplicates
    for entity in query:
        session.query(Entity).filter(Entity.pubmed_id == entity.pubmed_id, Entity.id != entity.id).delete(synchronize_session=False)

    # Commit the changes
    session.commit()

def get_info(keyword):
    pass