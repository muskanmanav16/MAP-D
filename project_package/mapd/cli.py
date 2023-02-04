import logging

import click
import uvicorn

from mapd.Database import Database, Utilapi

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
@click.group()
def main():
    """Entry method."""
    pass

@main.command()

def build_abstract_database():
    Database.build_database()
    Database.add_abstract_to_database()

@main.command():
def get_abstracts():

@main.command():
def entity_dict():
    Utilapi.get_abstracts()

@main.command()
@click.argument('pmid')
def get_abstract_info(pmid: str):
    """Retrieves identifier information for a given PubMed id."""

    entries_dict = Database.get_abstract_info(pubmed_id = pmid)

    return entries_dict
@main.command()
@click.option('-p', '--ppi', default=None, help="A CSV file containing PPIs.")
@click.option('-n', '--nodes', default=None, help="A TSV file containing defined nodes of a network.")
@click.option('-e', '--edges', default=None, help="A TSV file containing defined edges of a network.")
@click.option('-v', '--verbose', default=False, is_flag=True, help="Prints stats to STDOUT.")
@click.option('-r', '--enrich', default=False, is_flag=True, help="Enrich the graph with RNA and DNA molecules.")
@click.option('-o', '--output', default=None, help="File path to save summary stats.")
def stats(ppi: str, nodes: str, edges: str, output: str, verbose: bool, enrich: bool):
    """Generates node/edge lists from a PPI file."""
    sa = Statistics(node_list=nodes, edge_list=edges, ppi_file=ppi, enrich=enrich)
    sa.summary_statistics()

    if verbose:
        click.echo(sa.sum_stats)

    if output:
        sa.export_stats(output)
@main.command()
def serve():
    """Start the API web server."""
    uvicorn.run("project_package.frontend.Frontend_progress.run:app")


if __name__ == "__main__":
    main()
