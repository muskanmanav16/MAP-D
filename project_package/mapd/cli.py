import click
import uvicorn
import pandas as pd
from mapd.Database import Database
from mapd.NER import EntityPrediction
from mapd.utils import query_database
from mapd import DB_PATH
from pathlib import Path
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Instantiate Database class
db = Database()

@click.group()
def main():
    """Entry method"""
    pass


@main.command()
def build_db():
    """Builds database and populates it with abstract records using either cached/newly downloaded files."""
    db = Database()
    db.add_entity_data()
    click.echo("Built database of PubMed abstracts at {}.".format(DB_PATH))


@main.command()
def rebuild_db():
    """Rebuilds database from scratch after dropping any existing tables, and adds associated entities to Entity table"""

    db.rebuild_database()
    db.add_entity_data()
    click.echo("Rebuilt database of PubMed abstracts from scratch at {}.".format(DB_PATH))
    click.echo("Entity table of database populated with entities for each abstract.")


@main.command()
@click.option('-r', '--row_wise_results', default=False, is_flag=True, help="option to print dict entries row by row")
def get_entity_dict(row_wise_results: bool):
    """Returns a dictionary of entities in the table, with each entry containing
    an entity (key) and its label (value) """

    if db.raw_entity_data:
        entity_dict = db.raw_entity_data
    else:
        entity_dict = db.get_entity_dict()

    click.echo('Dictionary of entities (keys) and their labels (values)')
    if row_wise_results:
        for entry in entity_dict.items():
            click.echo(entry)
    else:
        click.echo(entity_dict)


@main.command()
@click.argument('pmid')
def get_abstract_info(pmid: int): # can test with PMID 36316711
    """Retrieves information about an abstract given its PMID (PubMed ID).

    Parameters
    ----------
    pmid: int
        PubMed ID of abstract
    """

    entries_dict = db.get_abstract_info(pubmed_id=pmid)
    click.echo(entries_dict)

@main.command()
@click.argument('text')
def predict_entities(text: str):
    """Given a text, predicts entities using Scispacy NER model.
    Parameters
    ----------
    text: str
        Input text for entity prediction
    """

    entity_predictor = EntityPrediction(db.session)
    entities = entity_predictor.predict_entities(text)
    click.echo('Entities predicted for input text:')
    click.echo(entities)


@main.command()
@click.argument('keyword')
@click.argument('filepath')
@click.option('-s', '--start_date', default=None, help='start date for time range of desired query result')
@click.option('-e', '--end_date', default=None, help='end date for time range of desired query result')
def query_db(keyword: str, filepath: str, start_date=None, end_date=None):
    # can test with: python cli.py query-db 'pancreatic' 'results2.csv' -s 2021-01-03 -e 2021-09-02
    """Queries database for keyword and (optionally) date range, and saves results to file at specified address.
    Parameters
    ----------
    keyword: str
        Keyword to query database with
    filepath: str
        File path for query results .csv file
    """

    results = query_database(keyword, start_date, end_date)
    for result in results:
        result['abstract_text'] = result['abstract_text'].replace('<mark>' + keyword + '</mark>', keyword)

    df = pd.DataFrame(results)
    df.to_csv(filepath, index=False)

    timerange = ' '
    if start_date and end_date:
        timerange = 'in time range ' + start_date + ' - ' + end_date
    click.echo('Results file for query "{}" {} stored at {}'.format(keyword, timerange, filepath))


@main.command()
def serve():
    """Starts web server."""

    uvicorn.run("frontend/Frontend_progress/run:app")


if __name__ == "__main__":
    main()