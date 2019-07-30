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
- `workflow`:
    - `lsf-queue`: name of the LSF queue to use
    - `dir`: directory to use for LSF job files. Must be accessible from the login **and** processing nodes.
    
## Usage

```bash
$ ipr-update-proteins config.json [OPTIONS]
$ ipr-update-pronto config.json [OPTIONS]
```
