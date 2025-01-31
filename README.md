# pyinterprod

A centralised Python implementation of InterPro production procedures.

## Getting started

### Requirements:

- Python 3.11+, with packages `oracledb`, `mysqlclient`, `psycopg3`, 
  and `mundone` ([link](https://github.com/matthiasblum/mundone))
- `GCC` with the `sqlite3.h` header

### Installation

```bash
pip install .
```

## Configuration

The `pyinterprod` package relies on three configuration files:

- `main.conf`: contains database connection strings, paths to files provided by/to UniProtKB, and various workflow parameters.
- `members.conf`: contains path to files used to update InterPro's member databases (e.g. files containing signatures, HMM files, etc.).
- `analyses.conf`: contains settings for the InterProScan match calculation (`ipr-calc`).
 
All files can be renamed. `main.conf` is passed as a command line argument, and the paths to `members.conf` and `analyses.conf` are defined in `main.conf`.  

### main.conf

> The expected format for database connection strings is `user/password@host:port/service`. 
> For Oracle databases, `user/password@service` may work as well, depending on `tnsnames.ora`.

- **oracle**
    - **ipro-interpro**: connection string for the `interpro` user in the InterPro database
    - **ipro-iprscan**: connection string for the `iprscan` user in the InterPro database
    - **ipro-uniparc**: connection string for the `uniparc` user in the InterPro database
    - **iscn-iprscan**: connection string for the `iprscan` user in the InterProScan database
    - **iscn-uniparc**: connection string for the `uniparc` user in the InterProScan database
    - **unpr-goapro**: connection string for the GOA database
    - **unpr-swpread**: connection string for the Swiss-Prot database
    - **unpr-uapro**: connection string for the UniParc production database 
    - **unpr-uaread**: connection string for the UniParc database
- **postgresql**:
    - **pronto**: connection string
- **uniprot**:
    - **version**: release number (e.g. `2019_08`)
    - **date**: date for the *public* release (e.g. `18-Sep-2019`)
    - **swiss-prot**: path to Swiss-Prot flat file
    - **trembl**: path to TrEMBL flat file
    - **unirule**: path to file listing InterPro entries and member database signatures used in UniRule
    - **xrefs**: path to directory where to export InterPro cross-references (generated for UniProt)
- **emails**:
    - **server**: outgoing server (format: `host:port`)
    - **sender**: sender's email address (e.g. user running the workflow)
    - **aa**: email address of the Automatic Annotation team
    - **aa_dev**: email address of the Automatic Annotation development team
    - **interpro**: email address of the InterPro team
    - **uniprot_db**: email address of the UniProt database team
    - **uniprot_db**: email address of the UniProt production team
    - **unirule**: email address of the UniRule team (curators from EMBL-EBI, SIB, and PIR)
    - **sib**: email address of the Swiss-Prot team
- **misc**:
    - **analyses**: path to the `analyses.conf` config file
    - **members**: path to the `members.conf` config file
    - **scheduler**: scheduler and queue (format: `scheduler:queue`, e.g. `lsf:production`)
    - **pronto_url**: URL of the Pronto curation application
    - **data_dir**: directory where to store staging files
    - **match_calc_dir**: directory where to run InterProScan match calculation
    - **temporary_dir**: directory for temporary files
    - **workflows_dir**: directory for workflows SQLite files, and jobs' input/output files

    
### members.conf

Each section corresponds to a member database (or a sequence feature database), e.g.

```
[profile]
signatures =
```

Supported properties are:

| Name         | Description                                                                                               |
|--------------|-----------------------------------------------------------------------------------------------------------|
| `signatures` | Path to the source of database signatures.                                                                |
| `hmm`        | Path to an HMM file, used for databases that employ HMMER3-based models. Required when running `ipr-hmm`. |
| `fasta`      | Path to sequences used by models, in the FASTA format.                                                    |
| `members`    | Path to file containing the clan-signature mapping.                                                       |
| `go-terms`   | Path to file or directory of GO annotations. PANTHER and NCBIfam only.                                    |
| `summary`    | Path to file of summary information. CDD only.                                                            |
| `seed`       | Path to file of SEED alignments. Pfam only.                                                               |
| `full`       | Path to file of full alignments. Pfam only.                                                               |
| `clans`      | Path to file of clan information. Pfam only.                                                              |
| `mapping`    | Path to file of model-signature mapping. CATH-Gene3D only.                                                |
| `classes`    | Path to file of information about classes. ELM only.                                                      |
| `instances`  | Path to file of information about instances. ELM only.                                                    |


### analyses.conf

The `DEFAULT` section defines the defaults values for the following properties:

- `job_cpu`: number of processes to request when submitting a job.
- `job_mem`: the maximum amount of memory a job should be allowed to use (in MB).
- `job_size`: the number of sequences to process in each job.
- `job_timeout`: the number of hours a job is allowed to run for before being killed. Any value lower than 1 disable the timeout.

The default values can be overridden. For instance, adding the following block under the `DEFAULT` section ensure that MobiDB-Lite jobs timeout after 48 hours and that PRINTS jobs are allocated 16GB of memory:

```
[mobidb-lite]
job_timeout = 48

[prints]
job_mem = 16384
```

## Usage

### Protein update

Update proteins and matches to the latest private UniProt release.

```bash
$ ipr-uniprot [OPTIONS] main.conf
```

The optional arguments are:

- `-t, --tasks`: list of tasks to run, by default all tasks are run (see [Tasks](#ipr-uniprot-tasks) for a description of available tasks)
- `--dry-run`: do not run tasks, only list those about to be run

<a name="ipr-uniprot-tasks"></a>

#### Tasks

<table>
<thead>
<tr>
    <th>Name</th>
    <th>Description</th>
    <th>Dependencies</th>
</tr>
</thead>
<tbody>
<tr>
    <td>update-uniparc</td>
    <td>Import UniParc cross-references</td>
    <td></td>
</tr>
<tr>
    <td>taxonomy</td>
    <td>Import the latest taxonomy data from UniProt</td>
    <td></td>
</tr>
<tr>
    <td>update-ipm-matches</td>
    <td>Update protein matches from ISPRO</td>
    <td></td>
</tr>
<tr>
    <td>update-ipm-sites</td>
    <td>Update protein site matches from ISPRO</td>
    <td></td>
</tr>
<tr>
    <td>update-proteins</td>
    <td>Import the new Swiss-Prot and TrEMBL proteins, and compare with the current ones</td>
    <td></td>
</tr>
<tr>
    <td>delete-proteins</td>
    <td>Delete obsolete proteins in all production tables</td>
    <td>update-proteins</td>
</tr>
<tr>
    <td>check-proteins</td>
    <td>Track UniParc sequences (UPI) associated to UniProt entries that need to be imported (e.g. new or updated sequence)</td>
    <td>delete-proteins, update-uniparc</td>
</tr>
<tr>
    <td>update-matches</td>
    <td>Update protein matches for new or updated sequences, run various checks, and track changes in protein counts for InterPro entries</td>
    <td>update-ipm-matches, check-proteins</td>
</tr>
<tr>
    <td>update-fmatches</td>
    <td>Update protein matches for sequence features (e.g. MobiDB-lite, Coils, etc.)</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>export-sib</td>
    <td>Export Oracle tables required by the Swiss-Prot team</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>report-changes</td>
    <td>Report recent integration changes to the UniRule team</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>aa-iprscan</td>
    <td>Build the AA_IPRSCAN table, required by the Automatic Annotation team</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>xref-condensed</td>
    <td>Build the XREF_CONDENSED table for the Automatic Annotation team (contains representations of protein matches for InterPro entries)</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>xref-summary</td>
    <td>Build the XREF_SUMMARY table for the Automatic Annotation team (contains protein matches for integrated member database signatures)</td>
    <td>report-changes</td>
</tr>
<tr>
    <td>export-xrefs</td>
    <td>Export text files containing protein matches for the UniProt database team</td>
    <td>xref-summary</td>
</tr>
<tr>
    <td>notify-interpro</td>
    <td>Notify the InterPro team that all tables required by the Automatic Annotation team are ready, so we can take a snapshot of our database</td>
    <td>update-fmatches, aa-iprscan, xref-condensed, xref-summary</td>
</tr>
<tr>
    <td>swissprot-de</td>
    <td>Export Swiss-Prot descriptions associated to member database signatures in the public release of UniProt (i.e. the release we are updating *from*)</td>
    <td></td>
</tr>
<tr>
    <td>unirule</td>
    <td>Update the list of signatures used by UniRule, so InterPro curators are warned if they attempt to unintegrated one of these signatures.</td>
    <td></td>
</tr>
<tr>
    <td>update-varsplic</td>
    <td>Update splice variant matches</td>
    <td>update-ipm-matches</td>
</tr>
<tr>
    <td>update-sites</td>
    <td>Update residue annotations</td>
    <td>update-ipm-sites, update-matches</td>
</tr>
<tr>
    <td><a href="#ipr-pronto-tasks">Pronto</a></td>
    <td>Update the Pronto PostgreSQL table</td>
    <td>taxonomy, update-fmatches, swissprot-de, unirule</td>
</tr>
<tr>
    <td>send-report</td>
    <td>Send reports to curators, and inform them that Pronto is ready</td>
    <td>Pronto tasks</td>
</tr>
</tbody>
</table>

### Member database update

Update models and protein matches for one or more member databases.

Before running the update, this command must be repeated for each member database. `-n` is the name of the database (case-insensitive), `-d` is the release date (of the member database), and `-v` is the release version.

```bash
$ ipr-pre-memdb main.conf -n DATABASE -d YYYY-MM-DD -v VERSION
```

Then, the actual update can be run:

```bash
$ ipr-memdb [OPTIONS] main.conf database [database ...]
```

The optional arguments are:

- `-t, --tasks`: list of tasks to run, by default all tasks are run (see [Tasks](#ipr-memdb-tasks) for a description of available tasks)
- `--dry-run`: do not run tasks, only list those about to be run

<a name="ipr-memdb-tasks"></a>

##### Tasks

<table>
<thead>
<tr>
    <th>Name</th>
    <th>Description</th>
    <th>Dependencies</th>
</tr>
</thead>
<tbody>
<tr>
    <td>update-ipm-matches</td>
    <td>Update protein matches from ISPRO</td>
    <td></td>
</tr>
<tr>
    <td>load-signatures</td>
    <td>Import member database signatures for the version to update to</td>
    <td></td>
</tr>
<tr>
    <td>track-changes</td>
    <td>Compare signatures between versions (e.g. name, description, matched proteins)</td>
    <td>load-signatures</td>
</tr>
<tr>
    <td>delete-obsoletes</td>
    <td>Remove signatures that are not in the latest version of the member database(s)</td>
    <td>track-changes</td>
</tr>
<tr>
    <td>update-signatures</td>
    <td>Update metadata for existing signatures, and add new signatures</td>
    <td>delete-obsoletes</td>
</tr>
<tr>
    <td>update-matches</td>
    <td>Update and check matches in production tables</td>
    <td>update-ipm-matches, update-signatures</td>
</tr>
<tr>
    <td>update-varsplic</td>
    <td>Update splice variant matches</td>
    <td>update-ipm-matches, update-signatures</td>
</tr>
<tr>
    <td>persist-pfam-a</td>
    <td>Parse Pfam-A files and store relevant information (only when updating Pfam)</td>
    <td>update-ipm-matches, update-signatures</td>
</tr>
<tr>
    <td>persist-pfam-c</td>
    <td>Parse Pfam-C to store clan information (only when updating Pfam)</td>
    <td>update-ipm-matches, update-signatures</td>
</tr>
<tr>
    <td>update-features</td>
    <td>Update sequence features for non-member databases (e.g. MobiDB-lite, COILS, etc.)</td>
    <td>update-ipm-matches</td>
</tr>
<tr>
    <td>update-fmatches</td>
    <td>Update matches for sequence features</td>
    <td>update-features</td>
</tr>
<tr>
    <td>update-ipm-sites</td>
    <td>Update protein site matches from ISPRO</td>
    <td></td>
</tr>
<tr>
    <td>update-sites</td>
    <td>Update residue annotations (if updating a member database with residue annotations)</td>
    <td>update-ipm-sites, update-matches</td>
</tr>
<tr>
    <td><a href="#ipr-pronto-tasks">Pronto</a></td>
    <td>Update the Pronto PostgreSQL tables</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>send-report</td>
    <td>Send reports to curators, and inform them that Pronto is ready</td>
    <td>Pronto tasks</td>
</tr>
</tbody>
</table>

### Pronto

```bash
$ ipr-pronto [OPTIONS] main.conf
```

The optional arguments are:

- `-t, --tasks`: list of tasks to run, by default all tasks are run (see [Tasks](#ipr-pronto-tasks) for a description of available tasks)
- `--dry-run`: do not run tasks, only list those about to be run

<a name="ipr-pronto-tasks"></a>

#### Tasks

<table>
<thead>
<tr>
    <th>Name</th>
    <th>Description</th>
    <th>Dependencies</th>
</tr>
</thead>
<tbody>
    <tr>
        <td>go-terms</td>
        <td>Import publications associated to protein annotations</td>
        <td></td>
    </tr>
    <tr>
        <td>go-constraints</td>
        <td>Import GO taxonomic constraints</td>
        <td></td>
    </tr>
    <tr>
        <td>proteins-similarities</td>
        <td>Import UniProt general annotations (comments) on sequence similarities</td>
        <td></td>
    </tr>
    <tr>
        <td>proteins-names</td>
        <td>Import UniProt sequence names</td>
        <td></td>
    </tr>
    <tr>
        <td>databases</td>
        <td>Import database information (e.g. version, release date)</td>
        <td></td>
    </tr>    
    <tr>
        <td>proteins</td>
        <td>Import general information on proteins (e.g. accession, length, species)</td>
        <td></td>
    </tr>
        <tr>
        <td>init-matches</td>
        <td>Create the match table (empty)</td>
        <td></td>
    </tr>
    <tr>
        <td>export-matches</td>
        <td>Export protein matches for member database signatures</td>
        <td>init-matches</td>
    </tr>
    <tr>
        <td>insert-matches</td>
        <td>Insert protein matches for member database signatures</td>
        <td>export-matches</td>
    </tr>
    <tr>
        <td>insert-fmatches</td>
        <td>Insert protein matches for sequence features (AntiFam, etc.)</td>
        <td>init-matches</td>
    </tr>
    <tr>
        <td>index-matches</td>
        <td>Index and cluster the match table</td>
        <td>insert-matches, insert-fmatches</td>
    </tr>
    <tr>
        <td>insert-signature2proteins</td>
        <td>Associate member database signatures with UniProt proteins, UniProt descriptions, taxonomic origins, and GO terms</td>
        <td>export-matches, proteins-names</td>
    </tr>
    <tr>
        <td>index-signature2proteins</td>
        <td>Index the signature2proteins table</td>
        <td>insert-signature2proteins</td>
    </tr>
    <tr>
        <td>signatures</td>
        <td>Import and compare member database signatures</td>
        <td>databases, export-matches</td>
    </tr>
    <tr>
        <td>taxonomy</td>
        <td>Import UniProt taxonomy</td>
        <td></td>
    </tr>
    <tr>
        <td>structures</td>
        <td>Import structural matches</td>
        <td></td>
    </tr>
</tbody>
</table>

### InterProScan match calculation

```bash
$ ipr-calc main.conf [COMMAND] [OPTIONS]
```

The available commands (and their optional arguments) are:

- `import`: import sequences from the UniParc Oracle database
  - `--top-up`: import new sequences only
- `clean`: delete obsolete data
  - `-a, --analyses`: IDs of analyses to clean (default: all)
- `search`: scan sequences using InterProScan
  - `-l, --list`: list active analyses and exit
  - `-a, --analyses`: IDs of analyses to run (default: all)
  - `--concurrent-jobs`: maximum number of concurrently running InterProScan jobs (default: 1000)
  - `--max-jobs`: maximum number of jobs to run per analysis before exiting (default: disabled)
  - `--max-retries`: number of times a failed job is resubmitted (default: disabled)
  - `--keep none|all|failed`: keep input/output files (default: none)

#### Examples

Import new UniParc sequences:

```bash
ipr-calc main.conf import --top-up
```

Process jobs for analysis `42` only, allow each job to run three times (i.e. restart twice), but keep all temporary files, regardless of the job success/failure:

```bash
ipr-calc main.conf search -a 42 --max-retries 2 --keep all
``` 

Run 10 jobs per analysis, and keep failed jobs to investigate:

```bash
ipr-calc main.conf search --max-retries 10 --keep failed
```

### Clans update

Update clans and run profile-profile alignments.

```bash
$ ipr-clans [OPTIONS] main.conf database [database ...]
```

The optional arguments are:

- `-t, --threads`: number of alignment workers
- `-T, --tempdir`: directory to use for temporary files

### HMMs update

Load HMMs in the database.

```bash
$ ipr-hmms main.conf database [database ...]
```
