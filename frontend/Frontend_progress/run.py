import csv
import os

from flask import Flask, render_template, request, make_response, send_file
import pandas as pd
# from dummy import get_info
from mapd.Database import Database

app = Flask(__name__, static_folder='static')



@app.route('/')
def home():
    """ renders the template search_form.html """
    #return render_template('search_form.html')
    return render_template('index.html')
    # return render_template('combi.html')

###############################

@app.route('/test')
def ss_form():
    """ renders the template search_form.html """
    #return render_template('search_form.html')
    return render_template('combi.html')

#######################################

# search_form.html contains:
# a form with a text input field for the keyword
# and two date inputs for the start and end date
# and a submit button that triggers the search_results()


@app.route('/results', methods=['POST'])
def search_results():
    """ retrieves the keyword, start and end dates from the form. Calls the get_info()
    with these parameters to retrieve relevant abstracts from the DB.
    Renders the template search_results.html with the results passed as a parameter"""
    keyword = request.form['keyword']
    start_date = request.form['start_date']
    end_date = request.form['end_date']
    db = Database()  # creating an instance of the Database class
    results = db.query_database(keyword, start_date, end_date)  # db_path might be a parameter in get_info() - modify acc.

    # paginate
    # can change this value by modifying the results_per_page variable.
    # The page variable is obtained from the query parameters in the request (e.g. /results?page=2).
    # The start and end variables are used to slice the results list so that only the desired chunk is returned.
    results_per_page = 2
    page = request.args.get("page", 1)
    start = (page - 1) * results_per_page
    end = start + results_per_page
    results = results[start:end]

    return render_template('combi.html', results=results)


# search_results.html
# displays the search results in a table format
# has the hyperlink to pubmed pg by mapping to the PMID
# includes a form with a submit button to download the results
# and a script that uses javascript to toggle the visibility of the abstracts


# search_form.html has a submit button to download file
@app.route('/download', methods=['POST'])
def download_results():
    """ converts the results to a dataframe to generate a CSV file, which
    the user can download."""
    keyword = request.form['keyword']
    start_date = request.form['start_date']
    end_date = request.form['end_date']
    db = Database()  # creating an instance of the Database class
    results = db.query_database(keyword, start_date, end_date)
    for result in results:
        result['id'] = f'=HYPERLINK("https://www.ncbi.nlm.nih.gov/pubmed/{result["PMID"]}","{result["PMID"]}")'
        result['abstract_text'] = result['abstract'].replace('<mark>'+keyword+'</mark>', keyword)
    df = pd.DataFrame(results)
    df.to_csv('results.csv', index=False)
    return send_file('results.csv', as_attachment=True)


# if __name__ == '__main__':
#     app.run(debug=True)

if __name__ == '__main__':
    # Check for the FLASK_PORT environment variable
    flask_port = os.getenv('FLASK_PORT')
    print(flask_port)
    # If the FLASK_PORT environment variable is set, use it as the port for the app
    if flask_port:
        app.run(debug=True, port=flask_port, host="0.0.0.0")
    else:
        # If the FLASK_PORT environment variable is not set, use the default port of 5000
        app.run(debug=True, port=5000, host="0.0.0.0")



