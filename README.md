# MAP-D: Next-Generation Search Engine for Academics and Researchers



## Description

MAP-D is an innovative search engine designed to make the process of finding relevant biomedical papers easier and more efficient. 

If you're a student, professor, or researcher, you'll find MAP-D a game-changer in your research process. Say goodbye to the frustration of traditional search engines and try MAP-D today!

### Features:
- Utilizes NER technology to categorize and index biomedical papers
- Delivers relevant search results from a large corpus
- Ideal for students, professors, and researchers alike

Try MAP-D now and take the first step towards a smoother and more efficient research experience!
## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Visuals

#### CSS Animations and Transitions

![CSS Animations and Transitions ↑](/uploads/b35e62291ed8c7845b7514265d8129f6/css_animations_and_transitions.webm)



## Installation
Installing scispacy requires two steps, installing the library and intalling the models. 

To install the library, run:
```bash
pip install scispacy
```

To install a model (see our full selection of available models below), run a command like the following:

```bash
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_sm-0.5.1.tar.gz
```

Note: We strongly recommend that you use an isolated Python environment (such as virtualenv or conda) to install scispacy.
Take a look below in the "Setting up a virtual environment" section if you need some help with this.
Additionally, scispacy uses modern features of Python and as such is only available for **Python 3.6 or greater**.


#### Setting up a virtual environment

[Conda](https://conda.io/) can be used set up a virtual environment with the
version of Python required for scispaCy.  If you already have a Python
environment you want to use, you can skip to the 'installing via pip' section.

1.  [Follow the installation instructions for Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html?highlight=conda#regular-installation).

2.  Create a Conda environment called "scispacy" with Python 3.9 (any version >= 3.6 should work):

    ```bash
    conda create -n scispacy python=3.9
    ```

3.  Activate the Conda environment. You will need to activate the Conda environment in each terminal in which you want to use scispaCy.

    ```bash
    source activate scispacy
    ```

Now you can install `scispacy` and one of the models using the steps above.


Once you have completed the above steps and downloaded one of the models below, you can load a scispaCy model as you would any other spaCy model. For example:
```python
import spacy
nlp = spacy.load("en_core_sci_sm")
doc = nlp("Alterations in the hypocretin receptor 2 and preprohypocretin genes produce narcolepsy in some animals.")
```

#### Note on upgrading
If you are upgrading `scispacy`, you will need to download the models again, to get the model versions compatible with the version of `scispacy` that you have. The link to the model that you download should contain the version number of `scispacy` that you have.

## Available Models

To install a model, click on the link below to download the model, and then run 

```python
pip install </path/to/download>
```

Alternatively, you can install directly from the URL by right-clicking on the link, selecting "Copy Link Address" and running 
```python
pip install CMD-V(to paste the copied URL)
```

| Model          | Description       | Install URL
|:---------------|:------------------|:----------|
| en_core_sci_sm | A full spaCy pipeline for biomedical data with a ~100k vocabulary. |[Download](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_sm-0.5.1.tar.gz)|
| en_core_sci_md |  A full spaCy pipeline for biomedical data with a ~360k vocabulary and 50k word vectors. |[Download](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_md-0.5.1.tar.gz)|
| en_core_sci_lg |  A full spaCy pipeline for biomedical data with a ~785k vocabulary and 600k word vectors. |[Download](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_lg-0.5.1.tar.gz)|
| en_core_sci_scibert |  A full spaCy pipeline for biomedical data with a ~785k vocabulary and `allenai/scibert-base` as the transformer model. You may want to [use a GPU](https://spacy.io/usage#gpu) with this model. |[Download](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_scibert-0.5.1.tar.gz)|
| en_ner_craft_md|  A spaCy NER model trained on the CRAFT corpus.|[Download](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_ner_craft_md-0.5.1.tar.gz)|
| en_ner_jnlpba_md | A spaCy NER model trained on the JNLPBA corpus.| [Download](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_ner_jnlpba_md-0.5.1.tar.gz)|
| en_ner_bc5cdr_md |  A spaCy NER model trained on the BC5CDR corpus. | [Download](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_ner_bc5cdr_md-0.5.1.tar.gz)|
| en_ner_bionlp13cg_md |  A spaCy NER model trained on the BIONLP13CG corpus. |[Download](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_ner_bionlp13cg_md-0.5.1.tar.gz)|


## CLI and Containerization

CLI command:

To build the database with two tables, Abstract and entity, run the following command:

```bash
python cli.py build-db
```

Containerization:

Please refer to our Dockerfile and run the following commands:

```bash
docker build . -t mapdpkg:latest
```

```bash
docker run --name groupprojectplab2 -p 5000:5000 -e DOCKER_CONTAINER=1 -d mapdpkg:latest
```

To stop the container:

```bash
docker stop groupprojectplab2
```
## Usage

![Quick_walkthrough-MAP-D ↑](/uploads/e03eaaef943bac975c2f046c564a17be/Quick_walkthrough-MAP-D.webm)

## Support
For any assistance, you can write to us at our official email address: mapd@gmx.net. <br>
-or- <br>
You can use the contact form on the website.

![image](/uploads/658761253ed260a6c88f4b323ee0e936/image.png)

We would love to hear from you!

## Roadmap



## Authors and acknowledgment
This project is contributed by:

**M**uskan Manav [s0mumana@uni-bonn.de]() ; [manavm0](https://gitlab.informatik.uni-bonn.de/manavm0)

**A**stha Anand [s0asanan@uni-bonn.de]() ; [ananda0](https://gitlab.informatik.uni-bonn.de/ananda0)

**P**arinishtha Bhalla [parinishtha.bhalla@gmail.com]() ; [bhallap0](https://gitlab.informatik.uni-bonn.de/bhallap0)

**D**eepika Pradeep [s0deprad@uni-bonn.de]() ; [pradeepd0](https://gitlab.informatik.uni-bonn.de/pradeepd0)

And special Thanks to our Mentor [Bruce Schultz](https://gitlab.informatik.uni-bonn.de/bschultz) for helping us out through out the project.


## Citing

We use **ScispaCy** model '**_en_ner_bionlp13cg_md_**' to create our search engine MAP-D which predict Named Entity from Pubmed Abstracts.

![](data/scispacy_logo.png)

Citing the paper which uses this project [ScispaCy: Fast and Robust Models for Biomedical Natural Language Processing](https://www.semanticscholar.org/paper/ScispaCy%3A-Fast-and-Robust-Models-for-Biomedical-Neumann-King/de28ec1d7bd38c8fc4e8ac59b6133800818b4e29). Additionally, indicating the version and model of ScispaCy used in the research.

```
@inproceedings{neumann-etal-2019-scispacy,
    title = "{S}cispa{C}y: {F}ast and {R}obust {M}odels for {B}iomedical {N}atural {L}anguage {P}rocessing",
    author = "Neumann, Mark  and
      King, Daniel  and
      Beltagy, Iz  and
      Ammar, Waleed",
    booktitle = "Proceedings of the 18th BioNLP Workshop and Shared Task",
    month = aug,
    year = "2019",
    address = "Florence, Italy",
    publisher = "Association for Computational Linguistics",
    url = "https://www.aclweb.org/anthology/W19-5034",
    doi = "10.18653/v1/W19-5034",
    pages = "319--327",
    eprint = {arXiv:1902.07669},
    abstract = "Despite recent advances in natural language processing, many statistical models for processing text perform extremely poorly under domain shift. Processing biomedical and clinical text is a critically important application area of natural language processing, for which there are few robust, practical, publicly available models. This paper describes scispaCy, a new Python library and models for practical biomedical/scientific text processing, which heavily leverages the spaCy library. We detail the performance of two packages of models released in scispaCy and demonstrate their robustness on several tasks and datasets. Models and code are available at https://allenai.github.io/scispacy/.",
}
```

ScispaCy is an open-source project developed by [the Allen Institute for Artificial Intelligence (AI2)](https://allenai.org/). AI2 is a non-profit institute with the mission to contribute to humanity through high-impact AI research and engineering.

