import sys
import argparse

def getMethodList(report):
    methods = []

    for line in report:
        method_ac = line.split('\t')[0]
        methods.append(method_ac)
        
    return methods

def generateReportDict(report):
    reportDict = {}
    for line in report:
        lines = []
        line = line.strip()
        method_ac = line.split("\t")[0]
        if method_ac in reportDict:
            lines = reportDict[method_ac]
            lines.append(line)
            reportDict[method_ac] = lines
        else:
            lines.append(line)
            reportDict[method_ac] = lines
    return reportDict
        
def filterReport(methods:list, reportDict:dict):
    methodSet = set(methods)
    reportSet = set(reportDict.keys())
    filtered_methods = reportSet.difference(methodSet)

    filtered_report = []
    for method in filtered_methods:
        filtered_report.append(reportDict[method])
    return filtered_report

def filter_match_count(report1File:str, report2File:str, outFileName:str):
    with open(report1File, 'r') as report1, open(report2File,'r') as report2, open(outFileName, 'w') as outFile:
        methods = []
        reportDict = {}
        filtered_report = []

        print('Start processing...')

        methods = getMethodList(report2)
        reportDict = generateReportDict(report1)
        # print(len(reportDict.keys())

        filtered_report = filterReport(methods, reportDict)
        print('Done processing. Writing out to a file...')

        outFile.write("\t".join(["Method", "Entry", "Previous value", "New value", "Change (%)\n"]))
        for line in filtered_report:
            new_line = line[0] + '\n'
            outFile.write(new_line)

if __name__ == "__main__":
    """ This script filter out overlaps in reports generated from member db update.
    Lines from report1 will be filtered out if they are found in report2. For example:
    - report1File is match_count.tsv
    - report2File is swissprotDEchange.tsv
    the new report1File will be a subset of match_count.tsv and swissprotDEchange.tsv will be intact.

    This script assumes that your files are tab delimited and the first column is method ac.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("report1File", help="report match_count file")
    parser.add_argument("report2File", help="report swiss_DE file")
    parser.add_argument("outFileName", help="Output file name")

    args = parser.parse_args()

    filter_match_count(args.report1File, args.report2File, args.outFileName)



