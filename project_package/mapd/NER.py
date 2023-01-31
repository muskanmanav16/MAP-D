import scispacy
import os
import sqlalchemy
import spacy
import logging
import tensorflow as tf
from sqlalchemy import create_engine, Column, Integer, String, Table, MetaData,ForeignKey, DATE,text,select, inspect,func
from sqlalchemy_utils import database_exists
from sqlalchemy.orm import declarative_base, sessionmaker
from typing import Optional, Union
import en_ner_bionlp13cg_md    #The model we are going to use
from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR
from mapd.models import Base,Abstract,Entity
from tqdm import tqdm


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

#check the device name and use GPU if it is available
device_name = tf.test.gpu_device_name()
if device_name != '/device:GPU:0':
    print('GPU device not found. Using CPU instead.')
    tf.config.set_visible_devices([], 'GPU')
else:
    print('Found GPU at: {}'.format(device_name))

#DB connection
Session = sessionmaker(bind=engine)
session = Session()
Base.metadata.create_all(engine)

# Load the sciSpacy model
model_name="en_ner_bionlp13cg_md"
nlp = spacy.load(model_name)


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


entity_predictor = EntityPrediction(session)

# Get all abstracts in the database from Abstract table
abstracts = session.query(Abstract).all()

# Loop through each abstract and predict entities
for abstract in tqdm(abstracts, desc="Predicting entities for abstracts"):
    entities = entity_predictor.predict_entities(abstract.abstract_text)
    entity_predictor.insert_entities(abstract.pubmed_id, entities)



# TO DO
#include above code in more structured form
#check for any duplicated entries in entity and abstarct table
#task3--make function to fetch the data for frontend
#include some ideas and functions for cli & unit_tests