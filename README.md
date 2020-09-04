# pyinterprod

A centralised Python/C implementation of InterPro production procedures.

## Getting started

### Requirements:

- Python 3.6+, with packages `cx_Oracle` and `mundone` ([link](https://github.com/matthiasblum/mundone))
- `GCC` with the `sqlite3.h` header

### Installation

```bash
python setup.py install
```

## Configuration

The `pyinterprod` package include command line tools, which expect a JSON config file. Properties are described bellow:

- `database`: connection information
    - `dsn`: Oracle data source name (format: _host:port/service_, or simply _service_)
    - `users`: connection strings (format: _user/password_)
- `release`: information on the UniProtKB release being updated
    - `version`: release number (e.g. `2019_08`)
    - `date`: date for the **public** release (e.g. `18-Sep-2019`)
- `paths`: files and directories used during the procedure
    - `flat-files`: Swiss-Prot and TrEMBL flat files
    - `results`: directory for any non-shared data, or large temporary files 
    - `unirule`: file listing InterPro entries and member database signatures used in UniRule
    - `xrefs`: InterPro cross-references, generated for UniProt
- `email-notifications`: boolean enabling or disabling email notifications to the InterPro or UniProt teams
- `workflow`:
    - `lsf-queue`: name of the LSF queue to use
    - `dir`: directory to use for LSF job files. Must be accessible from the login **and** processing nodes.
    
## Usage

### Protein update

```bash
$ ipr-uniprot CONFIG.JSON [OPTIONS]
```

The optional arguments are:

- `-t, --tasks`: list of tasks to run, by default all tasks are run (see [Tasks](#protein-update-tasks) for a description of available tasks)
- `--dry-run`: do not run tasks, only list those about to be run
- `--resume`: skip successfully completed tasks. By default all tasks are run, even if already completed in the past. 

<a name="protein-update-tasks"/>

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
    <td>load-proteins</td>
    <td>Import the new Swiss-Prot and TrEMBL proteins, and compare with the current ones</td>
    <td></td>
</tr>
<tr>
    <td>update-proteins</td>
    <td>Delete obsolete proteins in all production tables</td>
    <td>load-proteins</td>
</tr>
<tr>
    <td>update-uniparc</td>
    <td>Import UniParc cross-references</td>
    <td></td>
</tr>
<tr>
    <td>check-crc64</td>
    <td>Check that CRC64 checksums in UniProt entries and UniParc cross-references are the same</td>
    <td>update-proteins, update-uniparc</td>
</tr>
<tr>
    <td>proteins2scan</td>
    <td>Track UniParc sequences (UPI) associated to UniProt entries that need to be imported (e.g. new or updated sequence)</td>
    <td>check-crc64</td>
</tr>
<tr>
    <td>import-matches</td>
    <td>Import protein matches from ISPRO</td>
    <td></td>
</tr>
<tr>
    <td>import-sites</td>
    <td>Import residue annotations from ISPRO</td>
    <td></td>
</tr>
<tr>
    <td>update-variants</td>
    <td>Update splice variant matches</td>
    <td>update-uniparc, import-matches</td>
</tr>
<tr>
    <td>update-matches</td>
    <td>Update protein matches for new or updated sequences, run various checks, and track changes in protein counts for InterPro entries</td>
    <td>import-matches, proteins2scan</td>
</tr>
<tr>
    <td>update-feature-matches</td>
    <td>Update protein matches for sequence features (e.g. MobiDB-lite, Coils, etc.)</td>
    <td>import-matches, proteins2scan</td>
</tr>
<tr>
    <td>update-sites</td>
    <td>Update residue annotations for CDD, and SFLD</td>
    <td>import-sites, update-matches</td>
</tr>
<tr>
    <td>aa-iprscan</td>
    <td>Build the AA_IPRSCAN table, required by the Automatic Annotation team.</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>report-unintegrated</td>
    <td>Report recent integration changes to the Automatic Annotation and UniRule teams</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>xref-summary</td>
    <td>Build the XREF_SUMMARY table for the Automatic Annotation team (contains protein matches for integrated member database signatures)</td>
    <td>report-unintegrated</td>
</tr>
<tr>
    <td>xref-condensed</td>
    <td>Build the XREF_CONDENSED table for the Automatic Annotation team (contains representations of protein matches for InterPro entries)</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>alert-interpro</td>
    <td>Notify the InterPro team that all tables required by the Automatic Annotation team are ready, so we can take a snapshot of our database</td>
    <td>aa-iprscan, xref-summary, xref-condensed, update-feature-matches</td>
</tr>
<tr>
    <td>dump-xrefs</td>
    <td>Export text files containing protein matches</td>
    <td>xref-summary</td>
</tr>
<tr>
    <td>match-counts</td>
    <td>Update legacy materialized views</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>export-sib</td>
    <td>Export Oracle tables required by the Swiss-Prot team</td>
    <td>update-matches</td>
</tr>
<tr>
    <td>taxonomy</td>
    <td>Import the latest taxonomy data from UniProt</td>
    <td></td>
</tr>
<tr>
    <td>signatures-descriptions</td>
    <td>Copy the Swiss-Prot descriptions associated to member database signatures in the previous release of UniProt</td>
    <td></td>
</tr>
<tr>
    <td>pronto</td>
    <td>Refresh the database schema used by Pronto</td>
    <td>update-matches, update-feature-matches, taxonomy, signatures-descriptions</td>
</tr>
<tr>
    <td>report-curators</td>
    <td>Send reports to curators, and inform them that Pronto is ready</td>
    <td>pronto</td>
</tr>
<tr>
    <td>unirule</td>
    <td>Import the InterPro entries and member database signatures used in UniRule</td>
    <td>pronto</td>
</tr>
</tobdy>
</table>

### Pronto

```bash
$ ipr-pronto CONFIG.JSON [OPTIONS]
```

The optional arguments are:

- `-t, --tasks`: list of tasks to run, by default all tasks are run (see [Tasks](#pronto-tasks) for a description of available tasks)
- `-o, --output`: output report for curators (default: swiss_de_families.tsv in the current working directory) 

Note that there is not `--resume` argument. When using the `-t, --tasks` option, **only the passed tasks are run**, regardless of the status of their respective dependencies.

<a name="pronto-tasks"/>

#### Tasks

<table>
<thead>
<tr>
    <th>Name</th>
    <th>Description</th>
    <th>Source</th>
    <th>Dependencies</th>
    </tr>
</thead>
<tbody>
<tr>
    <td>annotations</td>
    <td>Import protein annotations</td>
    <td>Gene Ontology Annotation</td>
    <td></td>
</tr>
<tr>
    <td>publications</td>
    <td>Import publications associated to protein annotations</td>
    <td>Gene Ontology Annotation</td>
    <td></td>
</tr>
<tr>
    <td>terms</td>
    <td>Import GO terms</td>
    <td>Gene Ontology Annotation</td>
    <td></td>
</tr>
<tr>
    <td>databases</td>
    <td>Import database information (e.g. version, release date)</td>
    <td>InterPro</td>
    <td></td>
</tr>
<tr>
    <td>taxa</td>
    <td>Import taxonomy information</td>
    <td>InterPro</td>
    <td></td>
</tr>
<tr>
    <td>comments</td>
    <td>Import Swiss-Prot comments</td>
    <td>UniProt</td>
    <td></td>
</tr>
<tr>
    <td>descriptions</td>
    <td>Import protein descriptions</td>
    <td>UniProt</td>
    <td></td>
</tr>
<tr>
    <td>enzymes</td>
    <td>Import Enzyme Commission (EC) numbers</td>
    <td>UniProt</td>
    <td></td>
</tr>
<tr>
    <td>proteins</td>
    <td>Import Swiss-Prot and TrEMBL proteins</td>
    <td>InterPro</td>
    <td></td>
</tr>
<tr>
    <td>signatures</td>
    <td>Import member database signatures</td>
    <td>InterPro</td>
    <td></td>
</tr>
<tr>
    <td>matches</td>
    <td>Import protein matches, then count the number of proteins matched by each member database signature</td>
    <td>InterPro</td>
    <td>signatures</td>
</tr>
<tr>
    <td>signatures-proteins</td>
    <td>Associate member database signatures with UniProt proteins, UniProt descriptions, taxonomic origins, and GO terms</td>
    <td></td>
    <td>descriptions, signatures, taxa, terms</td>
</tr>
<tr>
    <td>index</td>
    <td>Create various database indexes</td>
    <td></td>
    <td>signatures-proteins</td>
</tr>
<tr>
    <td>compare</td>
    <td>Evaluate the similarity between member database signatures based on protein overlaps, common UniProt descriptions, common taxonomic origins, and common GO terms</td>
    <td></td>
    <td>signatures-proteins</td>
</tr>
<tr>
    <td>report</td>
    <td>Evaluate the gained and lost Swiss-Prot descriptions for InterPro entries between the previous and current UniProt releases</td>
    <td></td>
    <td>index</td>
</tr>
<tr>
    <td>copy</td>
    <td>Copy schema</td>
    <td></td>
    <td>compare, index</td>
</tr>
</tbody>
</table>

### Member database update

```bash
$ ipr-memberdb CONFIG_MEMBER.JSON [OPTIONS]
```

The optional arguments are:

- `-t, --tasks`: list of tasks to run, by default all tasks are run (see [Tasks](#memberdb-update-tasks) for a description of available tasks)
- `--dry-run`: do not run tasks, only list those about to be run
- `--resume`: skip successfully completed tasks. By default all tasks are run, even if already completed in the past. 

<a name="memberdb-update-tasks"/>

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
    <td>populate-method-stg</td>
    <td>Track UniParc sequences (UPI) associated to UniProt entries that need to be imported</td>
    <td></td>
</tr>
<tr>
    <td>generate-old-report</td>
    <td>Compare old and new member database methods and generate the newSigStatsReport file</td>
    <td>populate-method-stg</td>
</tr>
<tr>
    <td>update-method</td>
    <td>Update method table for the updated member database</td>
    <td>populate-method-stg</td>
</tr>
<tr>
    <td>update-db-version</td>
    <td>Update member database version</td>
    <td>update-method</td>
</tr>
<tr>
    <td>method2descriptions</td>
    <td>Update SwissProt descriptions for each method</td>
    <td>update-method</td>
</tr>
<tr>
    <td>proteins2scan</td>
    <td>Track UniParc sequences (UPI) associated to UniProt entries that need to be imported</td>
    <td></td>
</tr>
<tr>
    <td>update-iprscan2dbcode</td>
    <td>Update IDs for the member databases to update</td>
    <td></td>
</tr>
<tr>
    <td>create-match-tmp</td>
    <td>Select the methods and corresponding proteins to be updated and generate match_counts_new report for each member database updated</td>
    <td>update-iprscan2dbcode, proteins2scan</td>
</tr>
<tr>
    <td>update-match</td>
    <td>Exchange partition between new and old matches</td>
    <td>create-match-tmp</td>
</tr>
<tr>
    <td>update-site-match</td>
    <td>Update site matches (only for CDD and SFLD member databases)</td>
    <td>update-match</td>
</tr>
<tr>
    <td>drop-match-tmp</td>
    <td>Delete unused matches table</td>
    <td>update-site-match</td>
</tr>
<tr>
    <td>pronto</td>
    <td>Refresh the database schema used by Pronto</td>
    <td>update-match, method2descriptions</td>
</tr>
<tr>
    <td>pronto-copy</td>
    <td>Copy Pronto to another schema</td>
    <td>pronto</td>
</tr>
<tr>
    <td>report-curators</td>
    <td>Evaluate the gained and lost Swiss-Prot descriptions for InterPro entries between the previous and current member database releases and send an email containing an archive of the report files</td>
    <td>pronto</td>
</tr>

</tbody>
</table>
