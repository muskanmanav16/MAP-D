"""Collection of utility methods."""
# from mapd import DATA_DIR, engine, DB_PATH, PUBMED_DIR
# from Bio import Medline
# from tqdm import tqdm
# from sqlalchemy import create_engine, Column, Integer, String, Table, MetaData,ForeignKey, DATE,text,select, inspect,func
# from mapd.models import Base,Abstract,Entity
# from datetime import datetime
#
# from sqlalchemy.orm import sessionmaker
# Base.metadata.create_all(engine)
# Session = sessionmaker(bind=engine)
# session = Session()
# with engine.begin() as conn:
#     for file in tqdm(PUBMED_DIR.glob("*.xml"), desc="Processing Files"):
#         with open(file,encoding='utf-8') as handle:
#             records = Medline.parse(handle)
#             for record in records:
#                 pmid = record['PMID']
#                 abstract = record['AB']
#                 date_list = record.get('PHST', [''])
#                 if len(date_list) >= 4:
#                     date_p = date_list[3]
#                 else:
#                     date_p = ''
#                 try:
#                     print(date_p)
#                     date = datetime.strptime(date_p, '%Y/%m/%d %H:%M [pubmed]').date()
#                 except (KeyError, ValueError):
#                     date = None
#
#                 abstract_entry = Abstract(pubmed_id=pmid, abstract_text=abstract, date=date)
#                 session.add(abstract_entry)
#             session.commit()

