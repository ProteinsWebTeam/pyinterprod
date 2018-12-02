#include <Python.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sqlite3.h>


typedef struct entry_t {
    char accession[16];         // Accession ID
    char crc64[17];             // Sixty-four bit cyclic redundancy checksum
    short int is_reviewed;
    short int is_fragment;
    short int day;
    short int month;
    short int year;
    int taxon_id;               // Taxon ID
    int length;                 // Sequence length
    char identifier[17];        // Entry name
} entry_t;

void rtrim(char *str) {
    size_t n;
    n = strlen(str);
    while (n > 0 && isspace((unsigned char)str[n - 1])) {
        n--;
    }
    str[n] = '\0';
}

static PyObject *sprot_load(PyObject *self, PyObject *args) {
    char *src;
    char *dst;
    char *table;

    if (!PyArg_ParseTuple(args, "sss", &src, &dst, &table))
        return NULL;

    char zSql[1024] = "INSERT INTO ";
    strcat(zSql, table);
    strcat(zSql, " VALUES (?, ?, ?, ?, ?, ?, ?);");

    unsigned int num_entries = 0;

    FILE *fp = fopen(src, "r");
    if (fp == NULL)
        return NULL;

    sqlite3 *db;
    int rc = sqlite3_open(dst, &db);
    if (rc != SQLITE_OK) {
        fclose(fp);
        return NULL;
    }
    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, zSql, strlen(zSql), &stmt, NULL);
    if (rc != SQLITE_OK) {
        fclose(fp);
        sqlite3_close(db);
        return NULL;
    }

    char *errmsg;
    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errmsg);

    char buffer[1024];
    char month[4];
    char *str, *token, *saveptr, *ptr;
    char delimiters[] = " ";
    unsigned int i;
    entry_t e;

    while (fgets(buffer, 1024, fp)) {
        rtrim(buffer);

        if (strncmp(buffer, "ID", 2) == 0) {
            i = 0;
            for (str = buffer; ; str = NULL) {
                token = strtok_r(str, delimiters, &saveptr);
                if (token == NULL)
                    break;
                else if (i == 1)
                    strcpy(e.identifier, token);
                else if (i == 2) {
                    ptr = strstr(token, "Reviewed");
                    e.is_reviewed = ptr != NULL ? 1 : 0;

                } else if (i == 3)
                    e.length = atoi(token);
                i++;
            }
        } else if (strncmp(buffer, "AC", 2) == 0) {
            for (i = 0, str = buffer; ; str = NULL, i++) {
                token = strtok_r(str, delimiters, &saveptr);
                if (token == NULL)
                    break;
                else if (i && !strlen(e.accession)) {
                    // Remove the semi-colon
                    strncpy(e.accession, token, strlen(token)-1);
                    break;
                }
            }

        // } else if (strncmp(buffer, "DT", 2) == 0) {
        //     ptr = strstr(buffer, "sequence version");
        //     if (ptr != NULL) {
        //         e.day = atoi(&buffer[5]);
        //         strncpy(month, &buffer[8], 3);
        //         month[3] = 0;
        //
        //         if (strcmp(month, "JAN") == 0)
        //             e.month = 1;
        //         else if (strcmp(month, "FEB") == 0)
        //             e.month = 2;
        //         else if (strcmp(month, "MAR") == 0)
        //             e.month = 3;
        //         else if (strcmp(month, "APR") == 0)
        //             e.month = 4;
        //         else if (strcmp(month, "MAY") == 0)
        //             e.month = 5;
        //         else if (strcmp(month, "JUN") == 0)
        //             e.month = 6;
        //         else if (strcmp(month, "JUL") == 0)
        //             e.month = 7;
        //         else if (strcmp(month, "AUG") == 0)
        //             e.month = 8;
        //         else if (strcmp(month, "SEP") == 0)
        //             e.month = 9;
        //         else if (strcmp(month, "OCT") == 0)
        //             e.month = 10;
        //         else if (strcmp(month, "NOV") == 0)
        //             e.month = 11;
        //         else
        //             e.month = 12;
        //
        //         e.year = atoi(&buffer[12]);
        //     }
        } else if (strncmp(buffer, "DE   Flags:", 11) == 0) {
            ptr = strstr(buffer, "Fragment");
            if (ptr != NULL)
                e.is_fragment =  1;
        } else if (!e.is_fragment && strncmp(buffer, "FT   NON_TER", 12) == 0) {
            e.is_fragment = 1;
        } else if (strncmp(buffer, "OX", 2) == 0) {
            i = 0;
            for (str = buffer; ; str = NULL) {
                token = strtok_r(str, "=", &saveptr);
                if (token == NULL)
                    break;
                else if (i)
                    e.taxon_id = atoi(token);

                i++;
            }
        } else if (strncmp(buffer, "SQ", 2) == 0) {
            i = 0;
            for (str = buffer; ; str = NULL) {
                token = strtok_r(str, delimiters, &saveptr);
                if (token == NULL)
                    break;
                else if (i == 6)
                    strcpy(e.crc64, token);

                i++;
            }
        } else if (strncmp(buffer, "//", 2) == 0) {
            // test that rc (returned by sqlite_bind*) == SQLITE_OK
            sqlite3_bind_text(stmt, 1, e.identifier, -1, SQLITE_STATIC);
            sqlite3_bind_text(stmt, 2, e.accession, -1, SQLITE_STATIC);
            sqlite3_bind_int(stmt, 3, e.is_reviewed);
            sqlite3_bind_int(stmt, 4, e.is_fragment);
            sqlite3_bind_int(stmt, 5, e.length);
            sqlite3_bind_int(stmt, 6, e.taxon_id);
            sqlite3_bind_text(stmt, 7, e.crc64, -1, SQLITE_STATIC);
            sqlite3_step(stmt);
            sqlite3_reset(stmt);

            num_entries++;

            // Reset entry for next record
            memset(e.accession, 0, sizeof(e.accession));
            memset(e.identifier, 0, sizeof(e.identifier));
            memset(e.crc64, 0, sizeof(e.crc64));
            e.is_fragment = 0;
            e.is_reviewed = 0;
            e.taxon_id = 0;
            e.length = 0;
            // e.day = 0;
            // e.month = 0;
            // e.year = 0;
        }
    }

    fclose(fp);
    sqlite3_exec(db, "END TRANSACTION", NULL, NULL, &errmsg);
    sqlite3_finalize(stmt);
    sqlite3_close(db);
    return PyLong_FromUnsignedLong(num_entries);
}

static PyMethodDef SprotMethods[] = {
   {"load", sprot_load, METH_VARARGS,
    "Load an UniProtKB flat file into a SQLite database"},
   {NULL, NULL, 0, NULL }      /* Sentinel */
};

static struct PyModuleDef sprotmodule = {
    PyModuleDef_HEAD_INIT,
    "sprot",    /* name of module */
    NULL,       /* module documentation, may be NULL */
    -1,         /* size of per-interpreter state of the module,
                    or -1 if the module keeps state in global variables. */
    SprotMethods
};

PyMODINIT_FUNC PyInit_sprot(void) {
    return PyModule_Create(&sprotmodule);
}
