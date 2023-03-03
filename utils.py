import re
import numpy as np


def extract_level(tax, level_prefix, null_data="/"):
    g = re.findall(rf"({level_prefix}[^;|^$]*)(?=;|$)", tax)
    if len(g) == 0:
        return null_data

    return g[0].replace(level_prefix, "")


def backtrace_unassigned(row, unassigned_like_words, invalid_words):
    columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]
    try:
        valid_name = [val for val in row[columns].values if val not in invalid_words][
            -1
        ]
        row[columns] = [
            val if val not in invalid_words else f"Unclass. {valid_name}"
            for val in row[columns].values
        ]

        valid_name = [
            val for val in row[columns].values if val not in unassigned_like_words
        ][-1]
        row[columns] = [
            val if val not in unassigned_like_words else f"Unclass. {valid_name} {val}"
            for val in row[columns].values
        ]
    except:
        pass
    return row


def reassemble_taxon(row):
    columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]
    prefixes = ["d__", "p__", "c__", "o__", "f__", "g__"]
    return ";".join(["".join(pair) for pair in zip(prefixes, row[columns].values)])


def remove_duplicates(row):
    columns = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]
    try:
        values = [val for val in row[columns].values if "Unclass." not in val]
        unique_values, counts = np.unique(values, return_counts=True)
        duplicate_id = np.where(counts > 1)[0][0] if any(counts > 1) else None

        if duplicate_id is not None:
            count = counts[duplicate_id]
            duplicate_str = unique_values[duplicate_id]
            to_replace = count - 1
            new_vals = []
            for v in row[columns].values[::-1]:
                new_v = v
                if v == duplicate_str and to_replace > 0:
                    new_v = f"Unclass. {duplicate_str}"
                    to_replace -= 1
                new_vals.append(new_v)
            row[columns] = new_vals[::-1]
            
    except:
        pass
    return row


def backtrace_unassigned_fungi(row, unassigned_like_words, invalid_words):
    columns = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    try:
        valid_name = [val for val in row[columns].values if val not in invalid_words][
            -1
        ]
        row[columns] = [
            val if val not in invalid_words else f"Unclass. {valid_name}"
            for val in row[columns].values
        ]

        valid_name = [
            val for val in row[columns].values if val not in unassigned_like_words
        ][-1]
        row[columns] = [
            val if val not in unassigned_like_words else f"Unclass. {valid_name}_{val}"
            for val in row[columns].values
        ]
    except:
        pass
    return row


def reassemble_taxon_fungi(row):
    columns = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    prefixes = ["k__", "p__", "c__", "o__", "f__", "g__"]
    return ";".join(["".join(pair) for pair in zip(prefixes, row[columns].values)])
