from pathlib import Path
from typing import Optional, Union,List
from sqlalchemy import select, inspect
from sqlalchemy.orm import Session
from sqlalchemy_utils import database_exists
import pprint
from mapd import DATA_DIR, engine
from mapd.models import Base, Abstract,Entity
def get_abstract_info(pubmed_id: int) -> Optional[dict]:
    """Get abstract info  for a given pubmedid from the relational database."""
    session = Session(bind=engine)
    stmt = select(
        Abstract.pubmed_id,
        Abstract.Title,
        Abstract.date,
        Abstract.abstract_text,
        Entity.entity,
        Entity.labels
    ).filter_by(pubmed_id=pubmed_id).join(Entity, isouter=True)

    entries = session.execute(stmt).fetchall()
    ent={}
    for pubmed_id, title, date, abstract_text, entity, labels in entries:
        ent[labels]=entity
    for pubmed_id, title, date, abstract_text, entity, labels in entries:
        entries_dict = {
            "pubmed_id": pubmed_id,
            "Title": title,
            "date": date.strftime("%Y-%m-%d"),
            "abstract_text": abstract_text,
            "entities": ent
        }
    return entries_dict
pp = pprint.PrettyPrinter(indent=4) #just for printing
pp.pprint(get_abstract_info(36680181))
def get_info(keyword:str, start_date: str = None, end_date: str = None) -> List[dict]:

    pass


