import argparse
from BCBio import GFF

"""
Parse and compare two GFF files, reporting differences in features such as CDS, tRNA, ncRNA, regulatory regions, and rRNA.
"""

def parseGFF(gff_file):
    """
    Parse a GFF file and extract features ID.
    """
    gff = GFF.parse(open(gff_file))
    content = {}

    for rec in gff:
        for feature in rec.features:
            if feature.type not in content and feature.type in ["CDS", "tRNA", "ncRNA", "regulatory_region", "rRNA"]:
                content[feature.type] = []

            # Extract IDs of CDS (UniRef)
            if feature.type == "CDS":
                try:
                    uniref = [x for x in feature.qualifiers["Dbxref"] if x.startswith("UniRef:")]
                    content[feature.type].append(uniref[0])
                except:
                    pass

            # Extract IDs of tRNA (Sequence Ontology)
            elif feature.type == "tRNA":
                try:
                    so = [x for x in feature.qualifiers["Dbxref"] if x.startswith("SO:")]
                    content[feature.type].append(so[0])
                except:
                    pass

            # Extract IDs of ncRNA, regulatory regions, and rRNA (RFAM)
            elif feature.type in ["ncRNA", "regulatory_region", "rRNA"]:
                try:
                    rfam = [x for x in feature.qualifiers["Dbxref"] if x.startswith("RFAM:")]
                    content[feature.type].append(rfam[0])
                except:
                    pass

    return content

def summaryTable(gff1, gff2):

    # Prepare features data
    all_features = {}
    for type in set(gff1.keys()).union(set(gff2.keys())):
        all_features[type] = {"Common" : set(gff1[type]).intersection(set(gff2[type])) if type in gff1 and type in gff2 else [],
                              "Unique Gff1" : set(gff1[type]).difference(set(gff2[type])) if type in gff1 else [],
                              "Unique Gff2" : set(gff2[type]).difference(set(gff1[type])) if type in gff2 else []}
    
    # Write summary table
    print("Feature_Type\tCommon\tUnique_Gff1\tUnique_Gff2")
    for feature_type, data in all_features.items():
        print(f"{feature_type}\t{len(data['Common'])}\t{len(data['Unique Gff1'])}\t{len(data['Unique Gff2'])}")


if __name__ == "__main__":
    #Arguments parser
    parser=argparse.ArgumentParser(description="Compare two GFF files and report differences.")
    parser.add_argument("-1", dest="gff1", type=str, help="First GFF file", required=True)
    parser.add_argument("-2", dest="gff2",type=str, help="Second GFF file", required=True)
    args=parser.parse_args()

    # Main
    contentGFF1 = parseGFF(args.gff1)
    contentGFF2 = parseGFF(args.gff2)
    summaryTable(contentGFF1, contentGFF2)