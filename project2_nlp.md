# Programming Project 02 - Natural Language Processing

## Group Projects
These projects serve as a way for you to learn how to build a programmatic tool in a group setting. Often in this field, large packages and software are built in a collaborative effort, and knowing how to effectively communicate ideas, tasks, and workflow is an essential skill. Here, you will work as part of a group and demonstrate your ability to compose a cohesive (and working) tool comprised of components created by others *and* yourself.

## Deadline
**February 07, 2023 - 12:00/noon (UTC+1:00)**

## Submission Guidelines
1. Clone your repository to your local workstation/computer if you have not done so
2. Structure and work on your submission.
3. Commit your work to the git repository.
4. Create a git tag.
  * Tag name must be equivalent to "GroupProject".
  * Tag can be made via the command line or using the GitLab GUI
5. Be sure to include a PDF of your presentation in your repository once it is finished

## Package Requirements
* All code comprising the backend portion (i.e. the code responsible for downloading, parsing, and formatting the information) must be compiled as an installable Python package
* The package must contain the following:
    * A working CLI (see [CLI Requirements](#cli-requirements) below)
    * A clear and descriptive `README.md` that details what the package is, what it does, how it can be used, and examples of how to use the CLI
    * The necessary dependencies so that the package works immediately upon installation in a new virtual environment
    * Working unit tests that test at least 70% of the code in the package

## CLI Requirements
* Within your python package there must also be a working command line interface (CLI)
* CLI methods must contain proper documentation that is accessible via the command line
* CLI method documentation should contain:
    * Explanations of the arguments and options involved with the method
    * Brief description of what the method is used for

## Use of External Libraries
In general, one can make use of an external library or package that can aid in accomplishing a small subtask, such as a combinatorial problem, interface with an API, etc., but you cannot use a library or package capable of solving __all__ of your tasks. You are of course allowed to use modified code from your previous individual assignments (including that of PLAB1) where applicable. If you do choose to use an external resource to perform part of one of your tasks, it must be properly explained in the presentation. If you have any questions or concerns about whether a particular resource is allowed, please feel free to ask via email or issue.

## General Remarks
* The tasks are purposely written in such as manner as to require you, as a group, to figure out what tools are needed, what information needs to be gathered, and what resources should be used
* All code-based work is to be done in GitLab
* Use GitLab Issues to track and assign individual tasks and required work
* The software package (backend code) and web application (frontend) can be stored in separate folders in the root directory of your repository as shown here (you can rename these folders as you please):
```bash
├── frontend
└── project_package
```

## Grading (10 pts):

| Task | 1 | 2 | 3 | 4 | 5 |
|:---:|:---:|:---:|:---:|:---:|:---:|
| Points | 2 | 2 | 2 | 3 | 1 |

<div style="page-break-after: always"></div>

## Natural Language Processing - Introduction
Natural Language Processing (NLP) has been a huge field of research for several years, especially when it comes to language translations (e.g. deepl.com). Bioinformaticians began employing NLP for identifying interesting patterns in clinical texts and electronic health records in recent years, and several methods have been developed to help predict the development of certain diseases based on the data mined from these texts. A subfield of NLP is Named Entity Recognition (NER) in which individual terms are identified within text and tagged with specific labels. These can labels can be incredibly useful when it comes to searching through multiple texts since one doesn't need to use the exact term to find a relevant result. Most search engines today use some type of NER in their methods, and you will apply this technique to scientific papers to help researchers find articles of interest.


## Aims
The primary goal of this project is to develop a search engine using NER technology that will find relevant papers from a corpus.
1. __Create a corpus of scientific abstracts__ using the NCBI E-utilities API
2. __Apply NER technology__ to tag words within the abstracts with specific keywords
3. __Store the abstracts and their tagged metadata__ for faster searching
4. __Create a frontend__ which allows one to search the corpus using a keyword or phrase

## Tasks

### Task 1 - *Create a Corpus* (2 pts)
* Using tools available from NCBI, write code that downloads abstracts from PubMed from the last 3 years
    - Be sure to follow the guidelines of the API
    - Because we are only interested in biological publications, it may be helpful to filter the cached abstracts using their included metadata. There may be certain keywords or tags that can be used to select for relevant entries
    - If you find that the number of abstracts is too small to be useful, feel free to expand the number of abstracts collected from PubMed (e.g. from the last 4 years instead)

### Task 2 - *Tag the Abstracts using NER* (2 pts)
* Iterate through the abstracts and tag specific words or phrases using an NER model
    - While you are free to create your own model, it is recommended to find a pre-trained model or a library with a model included to use
* Create a way to store this information so that one can quickly access it
    - This can be done using cached files or a DB
    
### Task 3 - *Develop the Search Functionality* (2 pts)
* Generate functions that allows one to query your tagged corpus using a keyword or phrases and get informative results (see Task 3)
    - As an example, one should be able to input specific gene symbols (e.g. IL2, IFNA) or phrases (e.g. "cell cycle" or "apoptosis")\
* Write/modify methods for restricting the searches to a range of dates i.e. only search through abstracts within a certain range of time

### Task 4 - *GUI* (3 pts)
* Construct a web interface that allows one to search your corpus of abstracts using a specific keyword or phrase
    - This should include a box where users can input their keyword/phrase
    - There should be a feature included to filter the results to those abstracts within a specific range of dates
    - The returned results should include
        + The PMID of the abstract (with a hyperlink to its PubMed page)
        + The date of publication
        + A collapsable field with the entire abstract for users to read
        + A list of keywords/phrases in the abstract which were identified as being tagged with the user's search query
    - Include a button so that user's can download this information (you can choose the available format(s))

### Task 5 - *Containerize the Application* (1 pt)
* Create a `Dockerfile` in your project's root directory that successfully compiles your application into a Docker image
    * It should bundle your backend and frontend code together
    * You may choose any base image
    * All members of group should be in the image metadata using `LABEL`
    * The `ENTRYPOINT` should start your web application and the web app should be accessible using port mappings
