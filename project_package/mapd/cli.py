import click
import logging
import uvicorn
import pandas as pd
from mapd.Database import Database
from mapd.NER import EntityPrediction
from mapd.utils import query_database
from mapd import DB_PATH
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Instantiate Database class, add abstracts from cached/freshly downloaded files if they do not already exist:
db = Database()
db.add_abstract_to_database()
db.add_entity_data()

# Creating Click group:
@click.group()
def main():

    """Entry method"""
    pass

@main.command()
def rebuild_db():
    """Rebuilds database from scratch after dropping any existing tables, and adds associated entities to Entity table"""

    db.rebuild_database()
    db.add_abstract_to_database()
    db.add_entity_data()
    click.echo(
        """Rebuilt database of PubMed abstracts from scratch at {}.
        Entity table of database populated with entities for each abstract.""".format(DB_PATH))

@main.command()
@click.option('-r', '--row_wise_results', is_flag=True, default=False,
              help='option to print dictionary entries {entity:label} row by row')
def get_entity_dict(row_wise_results):
    """Returns a dictionary of entities in the table, with each entry containing an entity (key) and its label (value)"""

    entity_dict = db.get_entity_dict()
    click.echo('Dictionary of entities (keys) and their labels (values)')

    if row_wise_results:
        for entry in entity_dict.items():
            click.echo(entry)
    else:
        click.echo(entity_dict)

@main.command()
@click.argument('pmid', help ='PubMed ID of abstract for which information will be retrieved')
def get_abstract_info(pmid):
    """Retrieves information about an abstract given its PMID (PubMed ID)."""

    entries_dict = db.get_abstract_info(pubmed_id=pmid)
    click.echo(entries_dict)


@main.command()
@click.argument('abstract_id', help='ID of the associated abstract in the Abstract table of the database')
@click.argument('entities', help='Entities to be inserted into the database')
def insert_entities(abstract_id, entities):
    """Adds entities into the Entity table in the database for a given."""
    entity_predictor = EntityPrediction(db.session)
    entity_predictor.insert_entities(abstract_id, entities)
    click.echo('Entities inserted into Entity table for abstract {}'.format(abstract_id))

@main.command()
@click.argument('text', help='Input text for entity prediction')
def predict_entities(text):
    """Given a text, predicts entities using Scispacy NER model."""

    entity_predictor = EntityPrediction(db.session)
    entities = entity_predictor.predict_entities(text)
    click.echo('Entities predicted for input text:')
    click.echo(entities)

@main.command()
@click.argument('keyword')
@click.argument('filepath',help='File path for query results file')
@click.option('s','--start_date', default=None, help='start date for time range of desired query result')
@click.option('e','--end_date', default=None, help='end date for time range of desired query result')
def query_db(keyword: str, filepath, start_date=None, end_date=None):
    """Queries database for keyword and (optionally) date range, and saves results to file at specified address."""

    results = query_database(keyword, start_date, end_date)
    for result in results:
        result['id'] = f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/{result["id"]}","{result["id"]}")'
        result['abstract_text'] = result['abstract_text'].replace('<mark>' + keyword + '</mark>', keyword)

    df = pd.DataFrame(results)
    df.to_csv(filepath, index=False)

    timerange = ' '
    if start_date and end_date:
        timerange = 'in time range '+start_date+' and '+ end_date

    click.echo('Results file for query {}{} stored at {}'.format(keyword,timerange,filepath))


@main.command()
def serve():
    """Start the API web server."""
    uvicorn.run("group2.frontend.Frontend_progress.run:app") # this may need editing

if __name__ == "__main__":
    main()