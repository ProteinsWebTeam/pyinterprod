# pyinterprod

A centralised Python implementation of InterPro production procedures.

## Getting started

### Requirements:

- Python 3.7+, with packages `cx_Oracle`, `psycopg2`, and `mundone` ([link](https://github.com/matthiasblum/mundone))
- `GCC` with the `sqlite3.h` header

### Installation

```bash
python setup.py install
```

## Configuration

The `pyinterprod` package relies on two configuration files:

- `main.conf`: contains database connection strings, paths to files provided by/to UniProtKB, and various workflow parameters.
- `members.conf`: contains path to files used to update InterPro's member databases (e.g. files containing signatures, HMM files, etc.).
 
Both files can be renamed. `main.conf` is passed as an command line argument, and the path to `members.conf` is defined in `main.conf`.  

### main.conf

> The expected format for database connection strings is `user/password@host:port/service`. 
> For Oracle databases, `user/password@service` may work as well, depending on `tnsnames.ora`.

- **oracle**
    - **interpro**: connection string for the `interpro` user in the InterPro database
    - **iprscan**: connection string for the `iprscan` user in the InterPro database
    - **uniparc**: connection string for the `uniparc` user in the InterPro database
    - **goapro**: connection string for the GOA database
    - **swpread**: connection string for the Swiss-Prot database
    - **uaread**: connection string for the UniParc database
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
    - **members**: path to the `members.conf`
    - **pronto_url**: URL of the Pronto curation application
    - **data_dir**: directory where to store staging files
    - **lsf_queue**: name of the LSF queue to submit jobs to
    - **workflow_dir**: directory for temporary files (e.g. job input/output)
    
### members.conf

Each section corresponds to a member database (or a sequence feature database), e.g.

```
[pfam]
hmm =
members =
signatures =
```

The most common properties are:

- `hmm`: path the HMM file (only for databases using HMMER3-based models). Required when running `ipr-hmm`.
- `members`: path to the file containing the clan-signature mapping. Required when running `ipr-clans`.
- `signatures`: path to the source of the database signatures (database connection string for Pfam, file for other resources).

Additional properties:

- `summary`: CDD only, information on superfamilies
- `mapping`: Cath-Gene3D only, family-model mapping file


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
    <td>import-matches</td>
    <td>Import protein matches from ISPRO</td>
    <td>update-uniparc</td>
</tr>
<tr>
    <td>import-sites</td>
    <td>Import residue annotations from ISPRO</td>
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
    <td>import-matches, check-proteins</td>
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
    <td>import-matches</td>
</tr>
<tr>
    <td>update-sites</td>
    <td>Update residue annotations</td>
    <td>import-sites, update-matches</td>
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
    <td>import-matches</td>
    <td>Import protein matches from ISPRO</td>
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
    <td>import-matches, delete-obsoletes</td>
</tr>
<tr>
    <td>update-varsplic</td>
    <td>Update splice variant matches</td>
    <td>import-matches, delete-obsoletes</td>
</tr>
<tr>
    <td>update-features</td>
    <td>Update sequence features for non-member databases (e.g. MobiDB-lite, COILS, etc.)</td>
    <td>import-matches</td>
</tr>
<tr>
    <td>update-fmatches</td>
    <td>Update matches for sequence features</td>
    <td>import-features</td>
</tr>
<tr>
    <td>import-sites</td>
    <td>Import residue annotations from ISPRO (if updating a member database with residue annotations)</td>
    <td></td>
</tr>
<tr>
    <td>update-sites</td>
    <td>Update residue annotations (if updating a member database with residue annotations)</td>
    <td>import-matches, import-sites, update-matches</td>
</tr>
<tr>
    <td><a href="#ipr-pronto-tasks">Pronto</a></td>
    <td>Update the Pronto PostgreSQL table</td>
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
        <td>databases</td>
        <td>Import database information (e.g. version, release date)</td>
        <td></td>
    </tr>
    <tr>
        <td>annotations</td>
        <td>Import publications associated to protein annotations</td>
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
        <td>proteins</td>
        <td>Import general information on proteins (e.g. accession, length, species)</td>
        <td></td>
    </tr>
    <tr>
        <td>matches</td>
        <td>Import protein matches</td>
        <td>databases</td>
    </tr>
    <tr>
        <td>signature2proteins</td>
        <td>Associate member database signatures with UniProt proteins, UniProt descriptions, taxonomic origins, and GO terms</td>
        <td>proteins-names</td>
    </tr>
    <tr>
        <td>signatures</td>
        <td>Import member database signatures</td>
        <td>matches, signature2proteins</td>
    </tr>
    <tr>
        <td>taxonomy</td>
        <td>Import UniProt taxonomy</td>
        <td></td>
    </tr>
</tbody>
</table>

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
