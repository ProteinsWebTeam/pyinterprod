
import cx_Oracle
from .. import orautils, logger


def report_swissprot_changes(user, dsn, memberdb, prefix='swiss_de_report_'):
    con = cx_Oracle.connect(orautils.make_connect_string(user, dsn))
    cur = con.cursor()

    databases = {}
    analyses = []
    iprscanids = {'Q': 'HAMAP', 'D': 'SFLD', 'H': 'PFAM', 'P': 'PROSITE_PATTERNS', 'F': 'PRINTS', 'R': 'SMART', 'z': 'FUNFAM',
                  'N': 'TIGRFAM', 'M': 'PROSITE_PROFILES', 'Y': 'SUPERFAMILY', 'X': 'GENE3D', 'V': 'PANTHER', 'J': 'CDD', 'U': 'PIRSF'}

    for member in memberdb:
        cur.execute("""SELECT ID
                        FROM IPRSCAN.MV_SIGNATURE_LIBRARY_RELEASE
                        WHERE LIBRARY=:dbname and ROWNUM <= 2
                        ORDER BY ID DESC
                    """, dbname=iprscanids[member["dbcode"]])
        ids = [row[0] for row in cur]
        if ids:
            if len(ids) < 1:
                ids[1] = 0
        else:
            logger.info(
                f"member db {member['name']} not found in IPRSCAN.MV_SIGNATURE_LIBRARY_RELEASE")
            continue
        ids.sort()
        databases[member["dbcode"]] = (int(ids[0]), int(ids[1]))
        analyses += [int(ids[0]), int(ids[1])]

    query = """SELECT DISTINCT M.DBCODE, IPR.ANALYSIS_ID, IPR.METHOD_AC,
              D.TEXT, E.ENTRY_AC, E.NAME, E.ENTRY_TYPE
            FROM IPRSCAN.MV_IPRSCAN_MINI IPR
            INNER JOIN INTERPRO.PROTEIN_TO_SCAN PS
              ON IPR.UPI=PS.UPI
            INNER JOIN INTERPRO.PROTEIN P
              ON PS.PROTEIN_AC=P.PROTEIN_AC
            INNER JOIN INTERPRO_ANALYSIS_LOAD.PROTEIN_DESC PD
              ON P.PROTEIN_AC=PD.PROTEIN_AC
            INNER JOIN INTERPRO_ANALYSIS_LOAD.DESC_VALUE D
              ON PD.DESC_ID=D.DESC_ID
            INNER JOIN INTERPRO.METHOD M
              ON IPR.METHOD_AC=M.METHOD_AC
            INNER JOIN INTERPRO.ENTRY2METHOD EM
              ON M.METHOD_AC=EM.METHOD_AC
            INNER JOIN INTERPRO.ENTRY E
              ON EM.ENTRY_AC=E.ENTRY_AC
            WHERE IPR.ANALYSIS_ID IN({})
            AND P.DBCODE='S'
            AND P.FRAGMENT='N'
            """.format(','.join(str(item) for item in analyses))

    cur.execute(query)

    methods = {}
    for row in cur:
        dbcode = row[0]
        analysis_id = row[1]
        method_ac = row[2]
        descr = row[3]
        entry_ac = row[4]
        entry_name = row[5]
        entry_type = row[6]

        last_id, new_id = databases[dbcode]
        if dbcode in methods:
            db = methods[dbcode]
        else:
            db = methods[dbcode] = {}

        if method_ac in db:
            m = db[method_ac]
        else:
            m = db[method_ac] = {
                'acc': method_ac,
                'entry': (entry_ac, entry_name, entry_type),
                'analyses': {
                    last_id: set(),
                    new_id: set()
                }
            }

        m['analyses'][analysis_id].add(descr)

    cur.execute(
        """
        SELECT DBCODE, DBSHORT
        FROM INTERPRO.CV_DATABASE
        WHERE DBCODE IN({})
        """.format(','.join(f"'{item}'" for item in databases))
    )

    dbnames = dict(cur.fetchall())

    for dbcode in methods:
        last_id, new_id = databases[dbcode]
        lines = []
        for method_ac in methods[dbcode]:
            m = methods[dbcode][method_ac]

            entry_ac, entry_name, entry_type = m['entry']
            last_descrs = m['analyses'][last_id]
            new_descrs = m['analyses'][new_id]
            n_last = len(last_descrs)
            n_new = len(new_descrs)

            change = '{:.1f}'.format(n_new / n_last * 100) if n_last else 'N/A'

            gained = ' | '.join(new_descrs - last_descrs)
            lost = ' | '.join(last_descrs - new_descrs)

            # avoid writting in the file methods without gain or lost descriptions
            if gained or lost:
                lines.append((method_ac, entry_ac, entry_name,
                              entry_type, n_last, n_new, change, gained, lost))

        dbshort = dbnames[dbcode]

        with open(prefix + dbshort + '.tsv', 'wt') as fh:
            fh.write('Method\tEntry\tName\tType\t# of old descriptions\t# of new descriptions\tChange (%)\t'
                     'Descriptions gained\tDescriptions lost\n')

            for cols in sorted(lines, key=lambda x: (0 if x[3] == 'F' else 1, x[3], x[1])):
                fh.write('\t'.join(map(str, cols)) + '\n')
