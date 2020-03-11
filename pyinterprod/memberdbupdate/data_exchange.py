import time
import re
import cx_Oracle
import sys
import traceback
from .. import orautils


class Data_exchanger(object):
    def __init__(self, user, dsn, base_table):
        self.user, self.dsn = user, dsn

        if not base_table:
            print(
                "Error: base_table input param was empty, 'MATCH' or 'FEATURE_MATCH' is expected in production"
            )
            sys.exit(1)
        base_table = base_table.upper()
        if base_table == "MATCH" or base_table == "FEATURE_MATCH":
            # For InterPro production the base_table input parameter can only be:
            # (1) "MATCH":
            #     With MATCH_DBCODE_{} partition names, and MATCH_TMP and MATCH_NEW supporting tables.
            # (2) "FEATURE_MATCH":
            #     With FEATURE_MATCH_DBCODE_{} partition names, and FEATURE_MATCH_TMP and FEATURE_MATCH_NEW supporting tables
            self.base_table = base_table
        else:
            print(
                "base_table of 'MATCH' or 'FEATURE_MATCH' expected in production, but found '{0}' instead".format(
                    base_table
                )
            )
            sys.exit(1)
        self.connection = cx_Oracle.connect(
            orautils.make_connect_string(self.user, self.dsn)
        )
        self.cursor = self.connection.cursor()

    def check_dbcodes(self, dbcodes):
        # Confirm the DBCODEs look sensible
        for dbcode in dbcodes:
            query = 'select dbcode from {0} partition ("{0}_DBCODE_{1}") where rownum <2'.format(
                self.base_table, dbcode
            )
            self.cursor.execute(query)
            query_result = self.cursor.fetchone()[0]
            if query_result != dbcode:
                print(
                    "Error: Base table '{0}' and DB code '{1}' failed validation due to query '{2}'".format(
                        self.base_table, dbcode, query
                    )
                )
                print("Aborting!")
                sys.exit(1)

            if self.base_table == "MATCH" and dbcode in (
                "g",
                "j",
                "n",
                "q",
                "s",
                "v",
                "x",
            ):
                print(
                    "Error: Match base table '{0}' incompatible with (feature match) database code '{1}'".format(
                        self.base_table, dbcode
                    )
                )
                print("Aborting!")
                sys.exit(1)
            elif self.base_table == "FEATURE_MATCH" and dbcode not in (
                "g",
                "j",
                "n",
                "q",
                "s",
                "v",
                "x",
            ):
                print(
                    "Error: Feature match base table '{0}' incompatible with (match) database code '{1}'".format(
                        self.base_table, dbcode
                    )
                )
                print("Aborting!")
                sys.exit(1)

    def truncate_match_new(self):
        query = "TRUNCATE TABLE {0}_NEW".format(self.base_table)

        try:
            self.cursor.execute(query)
            print("{0}_NEW truncated.".format(self.base_table))
        except cx_Oracle.DatabaseError as exception:
            print("Failed to execute " + query)
            print(exception)
            exit(1)

    def exchange_partition(self, partitioned_tbl, dbcode):
        query = 'ALTER TABLE {0} EXCHANGE PARTITION("{1}_DBCODE_{2}") WITH TABLE {1}_NEW'.format(
            partitioned_tbl, self.base_table, dbcode
        )

        try:
            self.cursor.execute(query)
            print(
                "{0} exchanged with {1}_NEW.".format(partitioned_tbl, self.base_table)
            )
        except cx_Oracle.DatabaseError as exception:
            print("Failed to execute " + query)
            print(exception)
            exit(1)

    def finalise_table(self, table_name):
        self.rebuild_all_index(table_name)
        self.reenable_all_constraint(table_name)
        self.rebuild_partition(table_name)

    def rebuild_all_index(self, table_name):
        unusable_index = """SELECT index_name FROM user_indexes WHERE table_name = '{0}'
                            AND status = 'UNUSABLE'""".format(
            table_name
        )
        rebuild_index = "ALTER INDEX {0} REBUILD parallel 4"
        print("Re-building all indexes of " + table_name)
        self.execute_alter_stm(unusable_index, rebuild_index)

    def execute_alter_stm(self, select_stm, alter_stm):
        try:
            self.cursor.execute(select_stm)
            data = self.cursor.fetchall()
            if len(data) > 0:
                for row in data:
                    if len(row) == 2:
                        query = alter_stm.format(row[0], row[1])
                        self.cursor.execute(query)
                    else:
                        query = alter_stm.format(row[0])
                        self.cursor.execute(query)
            else:
                print("No unusable indexes or constraints.")
        except cx_Oracle.DatabaseError as exception:
            print("Failed to execute " + query)
            print(exception)
            exit(1)

    def reenable_all_constraint(self, table_name):
        disabled_constraint = """SELECT constraint_name FROM user_constraints WHERE table_name = '{0}'
                                 AND status = 'DISABLED'""".format(
            table_name
        )
        enable_constraint = "ALTER TABLE {0} ENABLE CONSTRAINT ".format(table_name)
        enable_constraint = enable_constraint + "{0}"
        print("Re-enabling all constraints of " + table_name)
        self.execute_alter_stm(disabled_constraint, enable_constraint)

    def rebuild_partition(self, table_name):
        unusable_partition = """SELECT p.index_name, partition_name FROM user_indexes i, user_ind_partitions p
                                WHERE  i.index_name = p.index_name and table_name = '{0}'
                                AND p.status = 'UNUSABLE'""".format(
            table_name
        )
        rebuild_partition = 'ALTER INDEX {0} REBUILD PARTITION "{1}" parallel 4'
        print("Re-building partitions of " + table_name)
        self.execute_alter_stm(unusable_partition, rebuild_partition)

    def get_count(self, dbcode):
        try:
            sqlLine = "SELECT COUNT(*) FROM {0} PARTITION({0}_DBCODE_{1})".format(
                self.base_table, dbcode
            )
            self.cursor.execute(sqlLine)
            for row in self.cursor:
                count = row[0]
            print(
                "\nTotal count at {0} table for member db {1} is: {2}".format(
                    self.base_table, dbcode, count
                )
            )
        except cx_Oracle.DatabaseError as exception:
            print("Failed to execute " + sqlLine)
            print(exception)
            exit(1)

    def close(self):
        self.cursor.close()
        self.connection.close()


def exchange_data(user: str, dsn: str, memberdb: list, base_table: str):

    data_exchanger = Data_exchanger(user, dsn, base_table)

    # Now perform the partition exchanges
    for member in memberdb:
        dbcode = member["dbcode"]
        print("Start exchanging data for member db: " + dbcode)
        print("---------------------------------------\n")

        # truncate match_new
        start_time = time.time()
        data_exchanger.truncate_match_new()
        print("Elapsed time in seconds: ")
        print(time.time() - start_time)
        print("\n")

        # rebuilding indexes and contraints of match_new
        start_time = time.time()
        data_exchanger.finalise_table("{0}_NEW".format(base_table))
        print("Elapsed time in seconds: ")
        print(time.time() - start_time)
        print("\n")

        # data exchange between match_tmp and match_new
        start_time = time.time()
        data_exchanger.exchange_partition("{0}_TMP".format(base_table), dbcode)
        print("Elapsed time in seconds: ")
        print(time.time() - start_time)
        print("\n")

        # rebuilding indexes and contraints of match_tmp
        start_time = time.time()
        data_exchanger.finalise_table("{0}_TMP".format(base_table))
        print("Elapsed time in seconds: ")
        print(time.time() - start_time)
        print("\n")

        # rebuilding indexes and contraints of match_new
        start_time = time.time()
        data_exchanger.finalise_table("{0}_NEW".format(base_table))
        print("Elapsed time in seconds: ")
        print(time.time() - start_time)
        print("\n")

        # data exchange between match and match_new
        start_time = time.time()
        data_exchanger.exchange_partition(base_table, dbcode)
        print("Elapsed time in seconds: ")
        print(time.time() - start_time)
        print("\n")

        # rebuilding indexes and contraints of match
        start_time = time.time()
        data_exchanger.finalise_table(base_table)
        print("Elapsed time in seconds: ")
        print(time.time() - start_time)
        print("\n")

        # match count
        start_time = time.time()
        data_exchanger.get_count(dbcode)
        print("Elapsed time in seconds: ")
        print(time.time() - start_time)
        print("\n")

    print("\n{0} table should have the new data now.".format(base_table))
    print(
        "Please compare the counts in this report to the one of {0}_TMP from CREATE_{0}_TMP.OUT to confirm this.".format(
            base_table
        )
    )

    data_exchanger.close()
