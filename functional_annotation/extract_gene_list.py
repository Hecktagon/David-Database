import csv

def extract_id_column(file, col_name):
    ids = set()
    sep = "\t" if file.endswith(".tsv") else ","

    with open(file, newline="") as f:
        reader = csv.reader(f, delimiter=sep)
        header = next(reader)
        idx = header.index(col_name)

        for row in reader:
            ids.add(row[idx])

    return ids


def ids_to_file(outfile, ids):
    with open(outfile, "w") as o:
        o.write("\n".join(ids) + "\n")


def write_gene_lists(infile, id_col, outfile):
    id_set = extract_id_column(infile, id_col)
    return id_set


if __name__ == "__main__":
    write_gene_lists("sample_rates_from_deuterater.csv", "analyte_id", "gene_list.txt")