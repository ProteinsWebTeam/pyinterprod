import time
import re
import cx_Oracle
import sys
import traceback
from .. import orautils,logger
from .match_tmp import add_indexes, add_constraints


class Data_exchanger(object):
    def __init__(self, user, dsn, base_table):
        self.user, self.dsn = user, dsn

        if not base_table:
            logger.info("Error: base_table input param was empty, 'MATCH' or 'FEATURE_MATCH' is expected in production")
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
            logger.info("base_table of 'MATCH' or 'FEATURE_MATCH' expected in production, but found '{0}' instead".format(base_table))
            sys.exit(1)
        self.connection = cx_Oracle.connect(orautils.make_connect_string(self.user, self.dsn))
        self.cursor = self.connection.cursor()

    def check_dbcodes(self, dbcodes:list):
        # Confirm the DBCODEs look sensible
        for dbcode in dbcodes:
            query = f"select dbcode from {self.base_table} partition ('{self.base_table}_DBCODE_{dbcode}') where rownum <2"
            self.cursor.execute(query)
            query_result = self.cursor.fetchone()[0]
            if query_result != dbcode:
                logger.info(f"Error: Base table '{self.base_table}' and DB code '{dbcode}' failed validation due to query '{query}'")
                logger.info("Aborting!")
                sys.exit(1)

            if self.base_table == "MATCH" and dbcode in ("g", "j", "n", "q", "s", "v", "x",):
                logger.info(f"Error: Match base table '{self.base_table}' incompatible with (feature match) database code '{dbcode}'")
                logger.info("Aborting!")
                sys.exit(1)
            elif self.base_table == "FEATURE_MATCH" and dbcode not in ("g", "j", "n", "q", "s", "v", "x",):
                logger.info(f"Error: Feature match base table '{self.base_table}' incompatible with (match) database code '{dbcode}'")
                logger.info("Aborting!")
                sys.exit(1)

    def recreate_match_new(self):
        table = f"{self.base_table}_NEW"
        orautils.drop_table(self.cursor, "INTERPRO", table , purge=True)
        
        query_create=f"CREATE TABLE {table} AS SELECT * FROM {self.base_table} WHERE 1=0"

        try:
            self.cursor.execute(query_create)
            logger.info(f"{table} created.")
        except cx_Oracle.DatabaseError as exception:
            logger.info(f"Failed to execute {query_create}")
            logger.info(exception)
            sys.exit(1)
        
        self.connection.commit()

        # add indexes
        try:
            add_indexes(self.cursor, self.connection, table, "MATCHN")
        except cx_Oracle.DatabaseError as e:
            logger.info(f"Failed to add indexes to {table} {e}")
            sys.exit(1)

        # add constraints
        try:
            add_constraints(self.cursor, self.connection, table, "MATCHN")
        except cx_Oracle.DatabaseError as e:
            logger.info(f"Failed to add constraints to {table} {e}")
            sys.exit(1)

    def truncate_match_new(self):
        query = f"TRUNCATE TABLE {self.base_table}_NEW"

        try:
            self.cursor.execute(query)
            logger.info(f"{self.base_table}_NEW truncated")
        except cx_Oracle.DatabaseError as exception:
            logger.info(f"Failed to execute {query}")
            logger.info(exception)
            sys.exit(1)
        

    def exchange_partition(self, partitioned_tbl:str, dbcode:str):
        query = f"ALTER TABLE {partitioned_tbl} EXCHANGE PARTITION({self.base_table}_DBCODE_{dbcode}) WITH TABLE {self.base_table}_NEW"

        try:
            self.cursor.execute(query)
            logger.info(f"{partitioned_tbl} exchanged with {self.base_table}_NEW.")
        except cx_Oracle.DatabaseError as exception:
            logger.info(f"Failed to execute {query}")
            logger.info(exception)
            exit(1)

    def finalise_table(self, table_name:str):
        self.rebuild_all_index(table_name)
        self.reenable_all_constraint(table_name)
        self.rebuild_partition(table_name)

    def rebuild_all_index(self, table_name:str):
        unusable_index = "SELECT index_name FROM user_indexes WHERE table_name = '{table_name}' AND status = 'UNUSABLE'"

        rebuild_index = "ALTER INDEX {0} REBUILD parallel 4"
        logger.info(f"Re-building all indexes of {table_name}")
        self.execute_alter_stm(unusable_index, rebuild_index)

    def execute_alter_stm(self, select_stm:str, alter_stm:str):
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
                logger.info("No unusable indexes or constraints.")
        except cx_Oracle.DatabaseError as exception:
            logger.info(f"Failed to execute {query}")
            logger.info(exception)
            exit(1)

    def reenable_all_constraint(self, table_name:str):
        disabled_constraint = f"SELECT constraint_name FROM user_constraints WHERE table_name = '{table_name}' AND status = 'DISABLED'"
        enable_constraint = f"ALTER TABLE {table_name} ENABLE CONSTRAINT "
        enable_constraint = enable_constraint + "{0}"
        logger.info("Re-enabling all constraints of {table_name}")
        self.execute_alter_stm(disabled_constraint, enable_constraint)

    def rebuild_partition(self, table_name:str):
        unusable_partition = f"SELECT p.index_name, partition_name FROM user_indexes i, user_ind_partitions p \
                                WHERE  i.index_name = p.index_name and table_name = '{table_name}' \
                                AND p.status = 'UNUSABLE'"

        rebuild_partition = 'ALTER INDEX {0} REBUILD PARTITION "{1}" parallel 4'
        logger.info("Re-building partitions of " + table_name)
        self.execute_alter_stm(unusable_partition, rebuild_partition)

    def get_count(self, dbcode:str):
        try:
            sqlLine = f"SELECT COUNT(*) FROM {self.base_table} PARTITION({self.base_table}_DBCODE_{dbcode})"
            self.cursor.execute(sqlLine)
            for row in self.cursor:
                count = row[0]
            logger.info(f"\nTotal count at {self.base_table} table for member db {dbcode} is: {count}")
        except cx_Oracle.DatabaseError as exception:
            logger.info(f"Failed to execute {sqlLine}")
            logger.info(exception)
            exit(1)

    def close(self):
        self.cursor.close()
        self.connection.close()


def exchange_data(user: str, dsn: str, memberdb: list, base_table: str):

    data_exchanger = Data_exchanger(user, dsn, base_table)

    # truncate match_new
    start_time = time.time()
    data_exchanger.recreate_match_new()
    logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

    # Now perform the partition exchanges
    for member in memberdb:
        dbcode = member["dbcode"]
        logger.info(f"Start exchanging data for member db: {dbcode}")

        # truncate match_new
        start_time = time.time()
        data_exchanger.truncate_match_new()
        logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

        # rebuilding indexes and contraints of match_new
        start_time = time.time()
        data_exchanger.finalise_table(f"{base_table}_NEW")
        logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

        # data exchange between match_tmp and match_new
        start_time = time.time()
        data_exchanger.exchange_partition(f"{base_table}_TMP", dbcode)
        logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

        # rebuilding indexes and contraints of match_tmp
        start_time = time.time()
        data_exchanger.finalise_table(f"{base_table}_TMP")
        logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

        # rebuilding indexes and contraints of match_new
        start_time = time.time()
        data_exchanger.finalise_table(f"{base_table}_NEW")
        logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

        # data exchange between match and match_new
        start_time = time.time()
        data_exchanger.exchange_partition(base_table, dbcode)
        logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

        # rebuilding indexes and contraints of match
        start_time = time.time()
        data_exchanger.finalise_table(base_table)
        logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

        # match count
        start_time = time.time()
        data_exchanger.get_count(dbcode)
        logger.info(f"Elapsed time in seconds: {time.time() - start_time}")

    logger.info(f"{base_table} table should have the new data now.")
    logger.info(f"Please compare the counts in this report to the one of {base_table}_TMP from CREATE_{base_table}_TMP.OUT to confirm this.")

    data_exchanger.close()
