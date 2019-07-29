# Open Source Economics Laboratory (OSE Lab) Boot Camp 2019

This public repository contains the training materials, tutorials, lecture notes, code, and problem sets for the seven-week Boot Camp of the Open Source Economics Laboratory (OSE Lab, [https://www.oselab.org/](https://www.oselab.org/)) at the University of Chicago, July 1 to August 9. The OSE Lab was founded by Dr. Richard W. Evans, Senior Lecturer and Associate Director at the University of Chicago M.A. Program in Computational Social Science. The OSE Lab is funded primarily from a 5-year grant from the Charles Koch Foundation. Part of this grant also included the creation of the Dynamic Analysis Center at the Baker Institute at Rice University, which is directed by John Diamond.

This year's OSE Lab Boot Camp has an exciting new partnership with the Econometric Society's Dynamic Structural Economics Summer School and Conference. July 8-9, 12-14 will be the [Dynamic Structural Economics Summer School](https://dseconf.org/dse2019course-program#) organized by [John Rust](https://editorialexpress.com/jrust/) (Georgetown University), [Bertel Schjerning](http://bschjerning.com/) (University of Copenhagen), and [Fedor Iskhakov](http://fedor.iskh.me/) (Australia National University). The Dynamic Structural Economics Conference, "[Applications in Industrial Organization, Marketing, and Business](https://editorialexpress.com/conference/DSE2019/program/DSE2019.html)" will take place on July 10-11 and is organized by [Sanjog Misra](https://www.chicagobooth.edu/faculty/directory/m/sanjog-misra) (Booth School of Business) and [Guenter Hitsch](https://www.chicagobooth.edu/faculty/directory/h/gunter-j-hitsch) (Booth School of Business). All OSE Lab students will be participants in both the DSE Summer School and in the conference.

This `README.md` serves as a syllabus and reference for the OSE Lab Boot Camp. This document has 12 sections.

1. [OSM Lab leadership](#1-ose-lab-leadership)
2. [Boot Camp schedule](#2-boot-camp-schedule)
3. [Instructions for installing the Anaconda distribution of Python](#3-instructions-for-installing-the-anaconda-distribution-of-python)
4. [Text editor suggestions](#4-text-editor-suggestions)
5. [PEP 8, docstring commenting, and module structure](#5-pep-8-docstring-commenting-and-module-structure)
6. [Using LaTeX](#6-using-latex)
7. [Git and GitHub.com tutorial](#7-git-and-gitHub-tutorial)
8. [Jupyter notebooks](#8-jupyter-notebooks)
9. [Python tutorials](#9-python-tutorials)
10. [Other Books](#10-other-books)
11. [C++ tutorials](#11-c-tutorials)
12. [References](#12-references)


## 1. OSE Lab Leadership

**Director: [Dr. Richard W. Evans](https://sites.google.com/site/rickecon/) [(rwevans@uchicago.edu)](mailto:rwevans@uchicago.edu).** The Director and Founder of the OSM Lab is [Richard W. Evans](https://sites.google.com/site/rickecon/). Dr. Evans is a Senior Lecturer and Associate Director of the M.A. Program in Computational Social Science at the University of Chicago. He is also a steering committee member of [QuantEcon](https://quantecon.org/).

**Director: [Dr. Simon Scheidegger](https://sites.google.com/site/simonscheidegger/home) [(simon.scheidegger@unil.ch)](mailto:simon.scheidegger@unil.ch).** This year's Boot Camp is Co-Organized by Dr. Simon Scheidegger. Simon is Assistant Professor in the Department of Finance at HEC Lausanne.

**Logistics and finances: Marty Evans, OSE Lab Administrative Director [(martye@uchicago.edu)](mailto:martye@uchicago.edu)**: Marty Evans is the Administrative Director of the OSE Lab. She coordinates the student and instructor travel, housing, and financial matters as well as other logistical needs.

**Instructors.** The OSE Lab has excellent instructors and presenters from economics and computation. Below is a list of this year's instructors and presenters in alphabetical order by last name for each of the three events in this year's boot camp. To see what they teach and when, check the respective [econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ), and [computation](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Computation) detailed schedule pages.

OSE Lab Senior Instructors and Presenters (in alphabetical order)
* [Scott Condie](https://scottcondie.github.io/), Brigham Young University
* [Richard Evans](https://sites.google.com/site/rickecon/), University of Chicago
* [Lars Hansen](http://larspeterhansen.org/), University of Chicago
* [Felix Kubler](https://www.bf.uzh.ch/en/persons/kuebler-felix), University of Zurich
* [Kerk Phillips](https://sites.google.com/site/kerkphillips/), Congressional Budget Office
* [John Rust](https://editorialexpress.com/jrust/), Georgetown University
* [Simon Scheidegger](https://sites.google.com/site/simonscheidegger/), HEC Lausanne
* [Anthony Smith, Jr.](http://www.econ.yale.edu/smith/), Yale University
* [Felix Tintelnot](http://felix-tintelnot.wikidot.com/), University of Chicago Economics
* [Thomas Winberry](http://faculty.chicagobooth.edu/thomas.winberry/), University of Chicago Booth School of Business

Graduate Instructors
* Marlon Azinovic, University of Zurich
* Rebekah Dix, Massachusetts Institute of Technology
* Jan Ertl, Oxford University

[DSE Summer School](https://dseconf.org/dse2019course-program#) Instructors (in alphabetical order)
* [Victor Aguirregabiria](http://individual.utoronto.ca/vaguirre/), University of Toronto
* [Øystein Daljord](https://faculty.chicagobooth.edu/oystein.daljord/), University of Chicago Booth School of Business
* [Paul Ellickson](http://paulellickson.com/), University of Rochester
* [Guenter Hitsch](https://www.chicagobooth.edu/faculty/directory/h/gunter-j-hitsch), University of Chicago Booth School of Business
* [Fedor Iskhakov](http://fedor.iskh.me/), Australia National University
* [Michael Keane](https://www.business.unsw.edu.au/our-people/mikekeane), University of New South Wales Business School
* [Robert A. Miller](https://www.cmu.edu/tepper/faculty-and-research/faculty-by-area/profiles/miller-robert.html), Carnegie Mellon University
* [Sanjog Misra](https://www.chicagobooth.edu/faculty/directory/m/sanjog-misra), University of Chicago Booth School of Business
* [Whitney Newey](https://economics.mit.edu/faculty/wnewey), MIT
* [Ariel Pakes](https://scholar.harvard.edu/pakes/home), Harvard University
* [John Rust](https://editorialexpress.com/jrust/), Georgetown University
* [Bertel Schjerning](http://bschjerning.com/), University of Copenhagen

[DSE Conference](https://editorialexpress.com/conference/DSE2019/program/DSE2019.html) Presenters and Organizers (in alphabetical order)
* [Bryan Bollinger](https://www.fuqua.duke.edu/faculty/bryan-bollinger), Duke University
* [Cheng Chou](https://chengchou.github.io/), University of Leicester
* [Giovanni Compiani](https://sites.google.com/berkeley.edu/giovannicompiani/home), University of California, Berkeley
* [Paul Ellickson](http://paulellickson.com/), University of Rochester
* [Guenter Hitsch](https://www.chicagobooth.edu/faculty/directory/h/gunter-j-hitsch), University of Chicago Booth School of Business
* [Jean Francois Houde](https://jfhoude.wiscweb.wisc.edu/), University of Wisconsin-Madison
* [Fedor Iskhakov](http://fedor.iskh.me/), Australia National University
* [Philipp Müller](https://www.business.uzh.ch/en/research/professorships/qba/members/mueller.html), University of Zurich
* [Sendhil Mullainathan](https://www.chicagobooth.edu/faculty/directory/m/sendhil-mullainathan), University of Chicago Booth School of Business
* [Charles Murray](https://murryecon.weebly.com/), Boston College
* [Harikesh Nair](https://www.gsb.stanford.edu/faculty-research/faculty/harikesh-s-nair), Stanford Graduate School of Business
* [Aviv Nevo](https://economics.sas.upenn.edu/people/aviv-nevo), University of Pennsylvania
* [Matthew Osborne](https://sites.google.com/site/matthewosborne/), University of Toronto
* [Harry Paarsch](https://business.ucf.edu/person/harry-paarsch/), University of Central Florida
* [Ariel Pakes](https://scholar.harvard.edu/pakes/home), Harvard University
* [Eduardo Souza-Rodrigues](http://individual.utoronto.ca/souza_rodrigues/), University of Toronto
* [John Rust](https://editorialexpress.com/jrust/), Georgetown University
* [Bertel Schjerning](http://bschjerning.com/), University of Copenhagen
* [Vira Semenova](https://economics.mit.edu/grad/vsemen), Harvard University
* [Timothy Simcoe](http://people.bu.edu/tsimcoe/), Boston University
* [K. Sudhir](https://som.yale.edu/faculty/k-sudhir), Yale School of Management
* [Ryan Webb](http://www.rotman.utoronto.ca/FacultyAndResearch/Faculty/FacultyBios/Webb), University of Toronto
* [Ruli Xiao](https://sites.google.com/site/iueconomicsrulixiao/), Indiana University


## 2. Boot Camp Schedule

The OSE Lab Boot Camp 2019 begins on Monday, July 1 and ends on Friday, August 9. Classes will be held Monday through Friday, 8am to noon in Saieh Hall, Room 247 at the University of Chicago. The exception is the Dynamic Structural Economics Workshop and Conference held during week 2, July 8-14. Workshop days on July 8-9 and 12-14 will be held at the University of Chicago Gleacher Center in Downtown Chicago. And the DSE conference will be held at the Booth School of Business. The curriculum consists of economic theory and computational methods with some math and statistics.

We have posted a one-page PDF of the [OSE Lab boot camp schedule](/OSE_Lab_topic_detail_2019.pdf) on this repository. Below is a summary schedule of topics. More detailed schedules are on the respective [econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ) and [computation](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Computation) detailed schedule pages.

### Week 1

| Date | Day | [Computation](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Computation#week-1) (8-10am) | [Econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ#week-1) (10am-noon) | [Computation](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Computation#week-1) (8am-noon) | Lunch speaker |
|:---:|:---:|:--- |:--- |:--- |:--- |
7-1  | M   | Numerical derivatives | Dynamic programming |              | Safety office |
7-2  | T   |              |                           | Python intro, standard library |     |
7-3  | W   | Numerical integration | Dynamic programming |              |     |
7-4  | Th  |              |                           |  |     |
7-5  | F   | Newton's method | Dynamic programming | Object oriented progr. |  |

### Week 2

Dynamic Structural Economics **Summer School** by the Econometric Society
* [Schedule](https://dseconf.org/dse2019course-program#), July 8, 9, 12, 13, 14

Dynamic Structural Economics **Conference** by the Econometric Society, "Applications of Structural Dynamic Models in Industrial Organization, Marketing, and Business Analytics"
* [Schedule](https://editorialexpress.com/conference/DSE2019/program/DSE2019.html), July 10, 11

### Week 3

| Date | Day | [Math](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Math#week-3) (8-10am) | [Econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ#week-3) (10am-noon) | [Computation](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Computation#week-3) (8am-noon) | Lunch speaker |
|:---:|:---:|:--- |:--- |:--- |:--- |
7-15  | M   | Measure theory | DSGE models |  | Supercomputer setup |
7-16  | T   |  |  | High performance computing |     |
7-17  | W   | Measure theory | DSGE models |       | [Jim Savage](https://schmidtfutures.com/person/jim-savage/), [Schmidt Futures](https://schmidtfutures.com/) |
7-18  | Th  |  |  | High performance computing |     |
7-19  | F   | Measure theory | DSGE models |  |  |

### Week 4

| Date | Day | [Econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ#week-4) (8-10am) | [Econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ#week-4) (10am-noon) | [Computation](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Computation#week-4) (8am-noon) | Lunch speaker |
|:---:|:---:|:--- |:--- |:--- |:--- |
7-22  | M   | Overlapping generations models | Firm Dynamics |  |  |
7-23  | T   |  |  | High performance computing | Marlon Azinovic, University of Zurich (slides, [paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3393482)) |
7-24  | W   | Overlapping generations models | Firm Dynamics |  |  |
7-25  | Th  |  |  | High performance computing | Lars Hansen, University of Chicago ([paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3008833)) |
7-26  | F   | Overlapping generations models | Firm Dynamics |  |  |

### Week 5

| Date | Day | [Econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ#week-5) (8-10am) | [Econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ#week-5) (10am-noon) | [Computation](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Computation#week-5) (8am-noon) | Lunch speaker |
|:---:|:---:|:--- |:--- |:--- |:--- |
7-29  | M   | Asset Pricing | Heterogeneous Agent Models |  |  |
7-30  | T   | Asset Pricing | Heterogeneous Agent Models |  |  |
7-31  | W   | Asset Pricing | Heterogeneous Agent Models |  |  |
8-1   | Th  |  |  | Pandas |  |
8-2   | F   |  |  | Conditioning, Stability, and Quasi-Newton |  |

### Week 6

| Date | Day | [Econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ#week-6) (8-10am) | [Econ](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Econ#week-6) (10am-noon) | [Computation](https://github.com/OpenSourceEcon/BootCamp2019/tree/master/Computation) (8am-noon) | Lunch speaker |
|:---:|:---:|:--- |:--- |:--- |:--- |
8-5  | M   | Machine Learning | International Models |  |  |
8-6  | T   | Machine Learning | International Models |  |  |
8-7  | W   | Machine Learning | International Models |  |  |
8-8  | Th  | Machine Learning |  |  |  |
8-9  | F   | Concluding discussion |  |  |  |

This will be an intensive six weeks. We expect that your attendance of lectures plus homework time might average 50-60 hours per week of work.

We have provided 6 areas of tutorials that you will benefit from reading and working through before the training. We will, of course, teach these things as we go through the material. But we will be able to proceed at a faster pace if the attendees are already familiar with most of the concepts below.

Pre-course Tutorial Areas

1. Instructions for installing the Anaconda distribution of Python
2. Text editor suggestions (Atom, Sublime Text 3, Vim)
3. PEP8, docstring commenting, and module structure
4. Git and GitHub tutorial
5. Jupyter Notebooks
6. Basic Python tutorials (data structures, logic, functions and modules, pandas, root finders and minimizers)


## 3. Instructions for installing the Anaconda distribution of Python

We will be using the [Python](https://www.python.org/) programming language and many of its powerful libraries for writing the code that will run most of the computational methods we will use during the Boot Camp. Using an open source language, such as Python, has the advantage of being free and accessible for anyone who wishes to learn these materials or contribute to these projects. Being open source also allows Python users to go into the source code of any function to modify it to suit one's needs.

We recommend that each participant download the Anaconda distribution of Python provided by [Anaconda, Inc.](https://www.anaconda.com/distribution/). We recommend the most recent stable version of Python, which is currently Python 3.7. This can be done from the [Anaconda download page](https://www.anaconda.com/distribution/) for Windows, Mac OSX, and Linux machines.


## 4. Text editor suggestions

In our recommended Python development workflow, you will write Python scripts and modules (`*.py` files) in a text editor. Then you will run those scripts from your terminal. You will want a capable text editor for developing your code. Many capable text editors exist, but we recommend three.

1. [Atom](https://atom.io)
2. [Sublime Text 3](https://www.sublimetext.com/)
3. [Vim](http://www.vim.org)

Atom and Vim are completely free. A trial version of Sublime Text 3 is available for free, but a licensed version is $80 (US dollars). In the following subsections, we give some of the details of each of the above three text editors.


### 4.1. Atom

[Atom](https://atom.io) is an open source text editor developed by people at GitHub.com. This editor has all the features of Sublime Text 3, but it also allows users full customizability. Further, it has been a while now that the users of Atom have surpassed the critical mass necessary to keep the editor progressing with the most cutting edge additions.

There are several packages you'll want to install with Atom.  Once Atom is installed, you can add packages by navigating Atom->Preferences->Install and then typing in the name of the package you would like to install.

For work with Python, we recommend the following packages be installed:

* MagicPython
* python-indent
* tabs-to-spaces
* minimap
* open-recent
* linter-python-pep8

For development with GitHub we recommend:

* merge-conflict

If using LaTex in this editor, the following packages are helpful:

* atom-latex
* latextools
* autocomplete-bibtex
* dictionary
* latexer
* pdf-view

In addition, if you are using a Mac, you will also want to download the [Skim](http://skim-app.sourceforge.net) PDF viewer to aid in displaying PDF files compiled from TeX with Atom.


### 4.2. Sublime Text 3

[Sublime Text 3](https://www.sublimetext.com) is the most widely used and versatile private software text editor. It has tremendous flexibility, as well as the polish of a piece of professional software. Sublime Text 3 will cost $80 for a license, although you can use a trial version indefinitely without charge while only having to suffer through frequent reminders to buy the full version.


### 4.3. Vim

[Vim](http://www.vim.org) is free and very powerful. Vim is the hard-core developer's text editor of choice. The learning curve for using vim is a little steeper than that of Atom and Sublime Text 3, but it also has some advantages for efficient programming. Vim has navigation that does not use a mouse or trackpad. Eventually, your fingers never leave your keyboard. Further, most terminals have Vim built in so you can more easily use Vim to edit scripts and modules on the fly with your terminal session.


## 5. PEP 8, docstring commenting, and module structure

Computer code executes some set of commands in an organized way. In every case, there are often many ways to execute a set of instructions--some ways more efficient than others. However, code has at least three functions.

1. Efficiently execute the task at hand.
2. Be accessible and usable to other programmers.
3. Be scalable and integrable with other projects and procedures.

Bill Gates is credited with the following plea for efficiency and parsimony in code writing.

> "Measuring programming progress by lines of code is like measuring aircraft building progress by weight."

Strong support for points (2) and (3) is Eagleson's Law.

> "Any code of your own that you haven't looked at for six or more months might as well have been written by someone else."

Because of the latter two characteristics, Python code has developed some conventions and best practices, some of which have been institutionalized in the [PEP 8--Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/) ("PEP" stands for Python Enhancement Proposals). Key examples PEP 8 Python coding conventions are the following.

* Indents should be 4 spaces (not tab)
* Limit all lines to a maximum of 79 characters long blocks of text being limited to 72 characters (Evans limits all his lines to 72 characters)
* Use a space after a comma
* Use a space before and after arithmetic operators

In the text editors Atom, Sublime Text 3, and Vim, you can install Linter packages that highlight areas of your code that break PEP 8 rules and tell you what the violation is.

There are fewer conventions in docstring structure, but we have developed some of our own that are outlined in the [PythonFuncs.ipynb](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonFuncs.ipynb) Jupyter notebook. See especially Sections 3 and 4 of the Jupyter notebook.


## 6. Using LaTeX

You will turn in all of your assignments by using the LaTeX document preparation platform. LaTeX produces documents with a sophisticated mathematical equation engine. Because LaTeX is standard in mathematical and theoretical document exposition, we will be using it in this class. The [LaTeX tutorial PDF chapter](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/LaTeX/LaTeX_tutorial.pdf) in the [`/Tutorials/LaTeX/`](https://github.com/OpenSourceMacro/BootCamp2019/tree/master/Tutorials/LaTeX) directory is a great reference for installing and running LaTeX. We have also included in that directory a template [`LaTeX_probset_template.tex`](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/LaTeX/LaTeX_probset_template.tex) as well as the PDf file [(`LaTeX_probset_template.pdf`)](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/LaTeX/LaTeX_probset_template.pdf) generated by compiling that `.tex` file.


## 7. Git and GitHub tutorial

We have included a tutorial on using [Git and GitHub.com](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/git_tutorial.pdf) in the [Tutorials](https://github.com/OpenSourceMacro/BootCamp2019/tree/master/Tutorials) directory of this repository. Git is a powerful version control software that comes natively installed on many machines and is widely used. GitHub.com is the most widely used online platform for hosting open source projects and integrating with Git software. Git has a significant learning curve, but it is essential for large collaborations that involve software development.

A more comprehensive Git resource is [*Pro Git*](https://git-scm.com/book/en/v2), by Chacon and Straub (2014). This book is open access, and is available online at [https://git-scm.com/book/en/v2](https://git-scm.com/book/en/v2). But Evans likes having it in his library in hard copy. This book is the difinitive guide with everything Git, and it has as its primary application the interaction between Git and GitHub. However, the workflow described in the tutorial above was hard to find in this Git book.


## 8. Jupyter Notebooks

[Jupyter notebooks](http://jupyter.org/) are files that end with the `*.ipynb` suffix. These notebooks are opened in a browser environment and are an open source web application that combines instructional text with live executable and modifyable code for many different programming platforms (e.g., Python, R, Julia). Jupyter notebooks are an ideal tool for teaching programming as they provide the code for a user to execute and they also provide the context and explanation for the code. We have provided a number of Jupyter notebooks in the [Tutorials](https://github.com/OpenSourceMacro/BootCamp2019/tree/master/Tutorials) folder of this repository.

These notebooks used to be Python-specific, and were therefore called iPython notebooks (hence the `*.ipynb` suffix). But Jupyter notebooks now support many programming languages, although the name still pays homage to Python with the vestigal "py" in "Jupyter". The notebooks execute code from the kernel of the specific programming language on your local machine.

Jupyter notebooks capability will be automatically installed with your download of the [Anaconda distribution](https://www.anaconda.com/download/) of Python. If you did not download the Anaconda distribution of Python, you can download Jupyter notebooks separately by following the instructions on the Jupyter [install page](http://jupyter.org/install).


### 8.1. Opening a Jupyter notebook

Once Jupyter is installed--whether through Anaconda or through the Jupyter website--you can open a Jupyter notebook by the following steps.

1. Navigate in your terminal to the folder in which the Jupyter notebook files reside. In the case of the Jupyter notebook tutorials in this repository, you would navigate to the `~/BootCamp2019/Tutorials/` directory.
2. Type `jupyter notebook` at the terminal prompt.
3. A Jupyter notebook session will open in your browser, showing the available `*.ipynb` files in that directory.
  *  In some cases, you might receive a prompt in the terminal telling you to paste a url into your browser.
4. Double click on the Jupyter notebook you would like to open.

It is worth noting that you can also simply navigate to the URL of the Jupyter notebook file in the GitHub repository on the web (e.g., [https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonReadIn.ipynb](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonReadIn.ipynb)). You can read the Jupyter notebook on GitHub.com, but you cannot execute any of the cells. You can only execute the cells in the Jupyter notebook when you follow the steps above and open the file from a Jupyter notebook session in your browser.


### 8.2. Using an open Jupyter notebook

Once you have opened a Jupyter notebook, you will find the notebook has two main types of cells: Markdown cells and Code cells. Markdown cells have formatted Jupyter notebook markdown text, and serve primarily to present context for the coding cells. A reference for the markdown options in Jupyter notebooks is found in the [Jupyter markdown documentation page](http://jupyter-notebook.readthedocs.io/en/latest/examples/Notebook/Working%20With%20Markdown%20Cells.html).

You can edit a Markdown cell in a Jupyter notebook by double clicking on the cell and then making your changes. Make sure the cell-type box in the middle of the top menu bar is set to `Markdown`. To implement your changes in the Markdown cell, type `Shift-Enter`.

A Code cell will have a `In [ ]:` immediately to the left of the cell for input. The code in that cell can be executed by typing `Shift-Enter`. For a Code cell, the  cell-type box in the middle of the top menu bar says `Code`.


### 8.3. Closing a Jupyter notebook

When you are done with a Jupyter notebook, you first save any changes that you want to remain with the notebook. Then you close the browser windows associated with that Jupyter notebook session. You should then close the local server instance that was opened to run the Jupyter notebook in your terminal window. On a Mac or Windows, this is done by going to your terminal window and typing `Cmd-C` or `Ctrl-C` and then selecting `y` for yes and hitting `Enter`.


## 9. Python tutorials

For this training, we have included in this repository six basic Python tutorials in the [`Tutorials`](https://github.com/OpenSourceMacro/BootCamp2019/tree/master/Tutorials) directory.

1. [PythonReadIn.ipynb](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonReadIn.ipynb). This Jupyter notebook provides instruction on basic Python I/O, reading data into Python, and saving data to disk.
2. [PythonNumpyPandas.ipynb](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonNumpyPandas.ipynb). This Jupyter notebook provides instruction on working with data using `NumPy` as well as Python's powerful data library `pandas`.
3. [PythonDescribe.ipynb](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonDescribe.ipynb). This Jupyter notebook provides instruction on describing, slicing, and manipulating data in Python.
4. [PythonFuncs.ipynb](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonFuncs.ipynb). This Jupyter notebook provides instruction on working with and writing Python functions.
5. [PythonVisualize.ipynb](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonVisualize.ipynb). This Jupyter notebook provides instruction on creating visualizations in Python.
6. [PythonRootMin.ipynb](https://github.com/OpenSourceMacro/BootCamp2019/blob/master/Tutorials/PythonRootMin.ipynb). This Jupyter notebook provides instruction on implementing univariate and multivariate root finders and unconstrained and constrained minimizers using functions in the [`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html) sub-library.

To further one's Python programming skills, a number of other great resources exist.

* The [official Python 3 tutorial site](https://docs.python.org/3/tutorial/)
* [QuantEcon.net](https://lectures.quantecon.org/py/) is a site run by [Thomas Sargent](http://www.tomsargent.com/) (NYU Stern) and [John Stachurski](http://johnstachurski.net/) (Australia National University). QuantEcon has a very large number of high-quality economics focused computational tutorials in Python. The first three sections provide a good introduction to Python programming.
* [Python computational labs](http://www.acme.byu.edu/2017-2018-materials/) of the Applied and Computational Mathematics Emphasis at Brigham Young University.
* [Code Academy's Python learning module](https://www.codecademy.com/learn/learn-python)

In addition, a number of excellent textbooks and reference manuals are very helpful and may be available in your local library. Or you may just want to have these in your own library. Lutz (2013) is a giant 1,500-page reference manual that has an expansive collection of materials targeted at beginners. Beazley (2009) is a more concise reference but is targeted at readers with some experience using Python. Despite its focus on a particular set of tools in the Python programming language, McKinney (2018) has a great introductory section that can serve as a good starting tutorial. Further, its focus on Python's data analysis capabilities is truly one of the important features of Python. Rounding out the list is Langtangen (2010). This book's focus on scientists and engineers makes it a unique reference for optimization, wrapping C and Fortran and other scientific computing topics using Python.


## 10. Other Books

Although the 2019 boot camp does not have extensive sections on pure mathematics anymore, we highly recommend the math book Humpherys, et al (2017). This book can be purchased through [SIAM](http://bookstore.siam.org/ot152/) or [Amazon](https://www.amazon.com/Foundations-Applied-Mathematics-Mathematical-Analysis/dp/1611974895).

The following books are recommendation that will not be required, but that are valuable in the library of a dynamic economist. We will also cover real analysis topics in measure theory. We recommend the following two books for their background on measure theory. Ok (2007) is a great real analysis book for economists. Kolmogorov and Fomin (1979) is a great real analysis book, the last four chapters of which are a good foundation for measure theory. Stokey and Lucas (1989, chap. 7) has an exposition of how measure theory gets used in dynamic stochastic macroeconomic theory. Stokey and Lucas (1989) is probably a good book for any macroeconomist to have in their library. On the economics side, some additional books that we like and recommend are Ljungqvist and Sargent (2012) and Adda and Cooper (2003).


## 11. C++ tutorials

Although we will be using Python for most of the Boot Camp, we will use some C++ for the computational labs taught by [Simon Scheidegger](https://sites.google.com/site/simonscheidegger/) on July 16, 18, 23, and 25. These computational labs treat high performance computing, parallel computing, and high dimensional approximation. Using supercomputing resources is much more flexible and accessible with C++.

[TODO: Include tutorial materials here.]


## 12. References

* Adda, Jerome and Russell Cooper, *Dynamic Economics: Quantitative Methods and Applications*, MIT Press (2003)
* Beazley, David M., *Python Essential Reference*, 4th edition, Addison-Wesley (2009).
* Chacon, Scott and Ben Straub, [*Pro Git: Everything You Need To Know About Git*](https://git-scm.com/book/en/v2), 2nd edition, Apress (2014).
* Humpherys, Jeffrey, Tyler J. Jarvis, and Emily J. Evans, *Foundations of Applied Mathematics, Volume 1: Mathematical Analysis*, SIAM (2017).
* Humpherys, Jeffrey, Tyler J. Jarvis, and Emily J. Evans, *Foundations of Applied Mathematics, Volume 2: Algorithm Design and Optimization*, SIAM (2018, forthcoming).
* Kolmogorov, A. N. and S.V. Fomin, *Introductory Real Analysis*, translated and edited by Richard A. Silverman, Dover Publications, Inc. (1970).
* Langtangen, Hans Petter, *Python Scripting for Computational Science*, Texts in Computational Science and Engineering, 3rd edition, Springer (2010).
* Ljungqvist, Lars and Thomas J. Sargent, *Recursive Methods in Macroeconomic Theory*, 3rd edition, MIT Press (2012).
* Lutz, Mark, *Learning Python*, 5th edition, O'Reilly Media, Inc. (2013)
* McKinney, Wes, *Python for Data Analysis*, 2nd edition, O'Reilly Media, Inc. (2018).
* Ok, Efe A., *Real Analysis with Economic Applications*, Princeton University Press (2007).
* Stokey, Nancy L., and Robert E. Lucas, Jr. with Edward C. Prescott, *Recursive Methods in Economic Dynamics*, Harvard University Press (1989).
