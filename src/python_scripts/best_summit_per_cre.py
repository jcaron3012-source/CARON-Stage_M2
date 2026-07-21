#!/usr/bin/env python3

'''Defines the best summit to use per CRE, with 3 methods: the summit nearest to the center, the best score or the median.
Returns the BED file with new columns defining the CRE summits.'''

import argparse
import pandas as pd
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cres",     default="/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/mus_links.bed",
                   help="CRE file (TSV, 1-based coordinates)")
    p.add_argument("--summits",  default="/Xnfs/lbmcdb/Semon_team/jcaron/files/peaks/by_cluster_by_molar/peaks_by_molar.bed",
                   help="Merged summits file (0-based, score in col5)")
    p.add_argument("--out",      default="/Xnfs/lbmcdb/Semon_team/jcaron/files/bedfiles/mus_with_best_summit.bed",
                   help="Output file")
    p.add_argument("--strategy", choices=["score", "center", "median"], default="score",
                   help="Selection criterion: highest score (default) or closest to CRE center")
    return p.parse_args()


def load_cres(path):
    df = pd.read_csv(path, sep="\t")
    # Rename seqnames to chr and convert 1-based start to 0-based
    df = df.rename(columns={"seqnames": "chr"})
    df["start0"] = df["start"] - 1
    df["end0"]   = df["end"]
    # Unique identifier per CRE
    df["cre_id"] = (df["chr"].astype(str) + ":"
                    + df["start0"].astype(str) + "-"
                    + df["end0"].astype(str))
    return df


def load_summits(path):
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["chr", "start", "end", "peak_id", "score"])
    # Extract cluster from peak_id, e.g. "mes-3-mx_peak_1" -> "mes-3-mx"
    df["cluster"] = df["peak_id"].str.replace(r"_peak_\d+$", "", regex=True)
    return df


def intersect(cres, summits):
    """Intersect CREs with summits by chromosome and coordinate overlap."""
    rows = []
    # Index summits by chromosome
    by_chr = {c: g for c, g in summits.groupby("chr")}

    for _, cre in cres.iterrows():
        chrom = str(cre["chr"])
        s, e = cre["start0"], cre["end0"]

        if chrom not in by_chr:
            continue

        hits = by_chr[chrom]
        # Summit is 1 bp: summit.start >= cre.start0 AND summit.end <= cre.end0
        mask = (hits["start"] >= s) & (hits["end"] <= e)
        matched = hits[mask]
        n_summits = int(mask.sum())

        for _, hit in matched.iterrows():
            rows.append({**cre.to_dict(), **{
                "n_summits": n_summits,
                "summit_chr": hit["chr"],
                "summit_start": hit["start"],
                "summit_end": hit["end"],
                "summit_id": hit["peak_id"],
                "summit_score": hit["score"],
                "cluster": hit["cluster"],
            }})

    return pd.DataFrame(rows)


def select_best(df, strategy, cres_original):
    '''Choose the best summit depending on the strategy'''
    
    if df.empty:
        print("WARNING: no intersections found.", file=sys.stderr)
        return df

    if strategy == "score":
        best = (df.sort_values("summit_score", ascending=False)
                  .drop_duplicates(subset="cre_id")
                  .reset_index(drop=True))
    elif strategy == "center":
        df = df.copy()
        df["cre_center"] = (df["start0"] + df["end0"]) // 2
        df["dist_center"] = abs(df["summit_start"] - df["cre_center"])
        best = (df.sort_values("dist_center")
                  .drop_duplicates(subset="cre_id")
                  .reset_index(drop=True))
    else:
        df = df.copy()
        # Compute the median summit position per CRE
        median_pos = (df.groupby("cre_id")["summit_start"]
                        .median()
                        .rename("median_summit_pos"))
        df = df.join(median_pos, on="cre_id")
        df["dist_median"] = abs(df["summit_start"] - df["median_summit_pos"])
        best = (df.sort_values("dist_median")
                  .drop_duplicates(subset="cre_id")
                  .reset_index(drop=True))


    # Append CREs with no summit (n_summits = 0)
    cres_with = set(best["cre_id"])
    cres_without = cres_original[~cres_original["cre_id"].isin(cres_with)].copy()
    cres_without["n_summits"] = 0
    cres_without["summit_chr"] = pd.NA
    cres_without["summit_start"] = pd.NA
    cres_without["summit_end"] = pd.NA
    cres_without["summit_id"] = pd.NA
    cres_without["summit_score"] = pd.NA
    cres_without["cluster"] = pd.NA

    result = pd.concat([best, cres_without], ignore_index=True)
    result = result.sort_values(["chr", "start0"]).reset_index(drop=True)
    return result


def main():
    args = parse_args()

    print(f"Loading CREs: {args.cres}")
    cres = load_cres(args.cres)
    print(f"  {len(cres):,} CREs loaded")

    print(f"Loading summits: {args.summits}")
    summits = load_summits(args.summits)
    print(f"{len(summits):,} summits loaded  |  {summits['cluster'].nunique()} clusters")

    print("Intersecting...")
    intersected = intersect(cres, summits)
    print(f"  {len(intersected):,} CRE x summit pairs")

    print(f"Selecting best summit (strategy: {args.strategy})")
    result = select_best(intersected, args.strategy, cres)

    n_with    = result["summit_id"].notna().sum()
    n_without = result["summit_id"].isna().sum()
    print(f"CREs with summit: {n_with:,}")
    print(f"CREs without summit: {n_without:,}")

    # Output columns: all original CRE columns + summit info
    keep_cre = [c for c in cres.columns if c not in ("start0", "end0", "cre_id")]
    out_cols  = keep_cre + ["n_summits", "summit_chr", "summit_start",
                             "summit_end", "summit_id", "summit_score", "cluster"]
    result[out_cols].to_csv(args.out, sep="\t", index=False)
    print(f"Output written: {args.out}")


if __name__ == "__main__":
    main()
