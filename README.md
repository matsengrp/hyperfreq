# hyperfreq

A Bayesian [APOBEC3G](http://en.wikipedia.org/wiki/APOBEC3G)-induced hypermutation analysis tool implemented in Python.


## CLI usage

The following examples illustrate some of the basic usage of hyperfreq.
Note that if one navigates to the root of this code base, these commands should all execute using the test data in `tests/data`.

### Analysis

To run an analysis, use the `hyperfreq analyze` command.
Running `hyperfreq analyze -h` will give you a full list of options.
Here are some examples to get you started and give you a rough sense of hyperfreq's capabilities.

    # Simple analysis comparing each sequence to a consensus sequence constructed from the entire alignment
    hyperfreq analyze tests/data/alignment.fasta

    # Instead, compare each sequence to the consensus for a cluster specified in a clusters file
    hyperfreq analyze tests/data/alignment.fasta -c tests/data/clusters.csv

    # Specify the reference sequence(s) you want each sequence to be compared to
    hyperfreq analyze tests/data/alignment.fasta -r tests/data/ref_seqs.fasta

These commands will all output data to a file named `hyperfreq_analysis.call.csv`.
The prefix and location of this file can be specified using the `-P` and `-o` flags.

### Multiple pattern analysis

By default, these analyses look for GG context hypermutation, suggestive of APOBEC3G activity.
One can specify multiple contexts for analysis using the `-p` or `--patterns` flag(s).
Pattern options include

 **pattern** | **associated with**
------------ | ------------------------------------------------------------------------
        `GG` | A3G activity
        `GA` | A3F (and other A3) activity in humans
        `GR` | combined A3G and A3F activity (as often observed in hypermutated HIV
        `GM` | rhesus macaque A3DE activity (as observed in XMRV and SFV infections) \*
        `GV` | combined rhesus A3DE and A3G activity

Note that R, M and V are IUPAC degenerate codes for A or G; A or C; and A, C or G, respectively.

When running multiple patterns, the `call.csv` file contains a column called `call_pattern` which represents the pattern in which the evidence of hypermutation appears to be strongest.
Other data in the `call.csv` file will contain counts, and statistics specifically for the pattern considered the call pattern.
If the `-F/--full-output` flag is specified, a separate file is output for each pattern analyzed (for example `hyperfreq_analysis.GG.csv`).

\* [Zhang et al](http://www.sciencedirect.com/science/article/pii/S0042682211004375)

   
### Splitting sequences for HM free alignments

Given an alignment, we can cut out sites/columns suspected of hypermutation by using the `split` command.
Running this command requires specifying an alignment and a CSV file with a column named `column`, specifying which positions in the alignment to be cut out.
In addition to the `hyperfreq_analysis.call.csv` file, `hyperfreq analyze` also produces a `hyperfreq_analysis.sites.csv` file which has such a column, as well as information regarding which sequences were hypermutated at which sequence positions.

For example,

    # for an alignment with hypermutated columns removed
    hyperfreq split alignment.fasta hypermutated_columns.csv

For more thorough usage, run `hyperfreq split -h` at the command line.


## Library usage

If you want to write your own scripts, you can do so by importing the appropriate modules

    from hyperfreq import Alignment
    from Bio import SeqIO

    # Create a hyperfreq alignment object
    seqs = SeqIO.parse('some_file.fasta', 'fasta')
    aln = Alignment(seqs)

    # Obtain an analysis generator which can be iterated over.
    analysis = aln.analyze()

    # Iterate over each sequence in the analysis, and do whatever you like!
    for seq_result in analysis:
        print seq_result['sequence'], "hm status:", seq_result['hm_pos']

It's also possible to define your own mutation patterns using the `MutPattern` and `MutPatternSet` classes.
It may be possible in the future to more flexibly specify patterns more flexibly via the CLI, but for now, doing so requires using this code base as a library in writing your own scripts.
If having this functionality available via the command line is important to you, please submit an [issue](https://github.com/fhcrc/hyperfreq/issues), and we'll see what we can do.


## Interactive usage

It is also now possible to run a `hyperfreq analyze` command with the `-N / --interactive` flag.
This causes the program to load up an interactive python session with the `analysis` generator, instead of simply writing the results to file.
You can enter `dir()` to see what namespaces and data have been included for you, and obviously load any libraries (numpy, biopython, etc.) that might be helpful.
Do note though that `analysis` is a generator, so if you need to make more than one pass through the results interactively, you'll need to throw them into an list.

    analysis = [result for result in analysis]


## Installing

Currently, Hyperfreq has only been tested on Linux systems, but it should be possible to get it set up on OSX fairly easily.
If you do get it set up on OSX or Windows and have any tips to share, please feel free to add a page to [the wiki](https://github.com/fhcrc/hyperfreq/wiki).

### Pre-requisites

Hyperfreq depends on the following python libraries

* `biopython`
* `fisher`
* [`betarat`](https://github.com/fhcrc/betarat)
* [`alnclst`](https://github.com/fhcrc/alnclst)

We recommend installing betarat and alnclst first, using the directions in the links above.
Since alnclst requires biopython, you should now only have to install fisher.
Assuming you have pip and the python-dev libraries installed (which you should after following the instructions above), you should be able to run

    # (sudo may not be necessary, depending on how you set up your python environment)
    sudo pip install fisher

### With that out of the way...

Now you should be able to download and install hyperfreq

    # I like to download things to a src directory in my home folder
    mkdir -p ~/src
    cd ~/src

    # Download and unzip the source code
    wget https://github.com/fhcrc/hyperfreq/archive/betarat-refactor.zip -O hyperfreq.zip
    unzip hyperfreq.zip
    cd hyperfreq

    # Install (sudo may not be necessary, depending on how you set up your python environment)
    sudo python setup.py install

And there you have it!
You can try running `hyperfreq -h` from the command line to test your installation.

### Problems?

If you have any trouble installing, please submit an [issue](https://github.com/fhcrc/hyperfreq/issues), so we can try to help and update the documentation.


