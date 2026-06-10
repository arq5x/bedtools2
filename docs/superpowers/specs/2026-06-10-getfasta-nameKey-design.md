# getfasta `-nameKey`: name FASTA headers from a GFF3 attribute

**Issue:** https://github.com/arq5x/bedtools2/issues/182
**Date:** 2026-06-10
**Status:** Design approved, pending implementation plan

## Problem

For GFF input, `parseGffLine` (`src/utils/bedFile/bedFile.h`) sets `bed.name = fields[2]`,
which is the GFF **type** column (`gene`, `mRNA`, `exon`, ...), not a meaningful name. The
informative identifiers live in the column-9 attributes string (`ID=...;Name=...`), stored in
`bed.fields[8]` but never consulted when building a FASTA header.

Consequently `getfasta -name` on a GFF emits useless headers like `>mRNA::chr1:99-200`. Issue
#182 asks for `getfasta` to instead pull `ID`/`Name` from the GFF3 attributes column.

## Constraint

The change must not alter the output of any existing command. Redefining `-name` / `-nameOnly`
for GFF would silently change output for current users, so the solution is a **new, opt-in flag**.

## Solution overview

Add a new flag to `getfasta`:

```
-nameKey <attr>   Use the value of GFF column-9 attribute <attr> for the FASTA header.
```

- Strictly additive. `-name`, `-name+`, `-nameOnly`, and the default behavior are untouched.
- Default output is **name-only**: `getfasta -bed a.gff -nameKey Name` -> `>geneA`.
- GFF3 `key=value;` attribute syntax only.

### Worked example

Input GFF3 record:

```
chr1  src  mRNA  100  200  .  +  .  ID=g1;Name=geneA
```

```
$ getfasta -fi g.fa -bed a.gff -nameKey Name
>geneA
ACGT...
```

## Detailed behavior

### CLI parsing (`src/fastaFromBed/fastaFromBedMain.cpp`)

- New variables: `bool useNameKey = false;`, `string nameKey;`.
- New branch consuming one argument (like `-fi`):

  ```cpp
  else if(PARAMETER_CHECK("-nameKey", 8, parameterLength)) {
      if ((i+1) < argc) {
          useNameKey = true;
          nameKey = argv[i + 1];
          i++;
      }
  }
  ```

- No prefix collision: the existing `PARAMETER_CHECK` macro requires an exact length match
  (`actualLen == paramLen`), so `-nameKey` (len 8) cannot match `-name` (5), `-name+` (6),
  or `-nameOnly` (9), and vice-versa.

### Mutual exclusivity (validation)

`-nameKey` is mutually exclusive with `-name`, `-name+`, and `-nameOnly`. If more than one of
these naming flags is supplied, print an error and exit before doing any work:

```
*****ERROR: -nameKey cannot be combined with -name, -name+, or -nameOnly. Choose one. *****
```

(Checked in `fastafrombed_main` after the parse loop, alongside the existing `haveFastaDb` /
`haveBed` checks.)

### Attribute extraction (`src/fastaFromBed/fastaFromBed.cpp`)

New helper, e.g. `bool getGffAttribute(const BED &bed, const string &key, string &value)`:

1. Obtain the attributes string: `bed.fields[8]` when `bed.fields.size() >= 9`; otherwise the
   record has no attributes column -> return `false`.
2. Split the attributes string on `;`.
3. For each token: trim leading/trailing whitespace, split on the first `=`. Compare the
   left side to `key` (case-sensitive). On match, set `value` to the right side, stripping a
   single pair of surrounding double-quotes if present, and return `true`.
4. If no token matches, return `false`.

Limitations (documented): GFF3 only (no GTF `key "value";`); no URL/percent-decoding of GFF3
escaped characters.

### Header construction (`ReportSeq` in `src/fastaFromBed/fastaFromBed.cpp`)

When `_useNameKey` is set:

- If `getGffAttribute` returns `true`, the name source is the extracted value (name-only).
- If it returns `false`, emit a per-record warning to stderr and fall back to the current
  `bed.name` (name-only):

  ```
  WARNING: attribute 'Name' not found for feature at chr1:99-200; using feature name.
  ```

- The existing `-s` strand-suffix path is unchanged, so `-nameKey Name -s` -> `>geneA(+)`.
- The run never aborts mid-stream because of a missing attribute.

The `Bed2Fa` constructor gains two parameters (`bool useNameKey`, `const string &nameKey`)
stored as `_useNameKey` / `_nameKey`.

## Files touched

- `src/fastaFromBed/fastaFromBedMain.cpp` — parse `-nameKey`, mutual-exclusivity check, help text.
- `src/fastaFromBed/fastaFromBed.h` — constructor signature, `_useNameKey` / `_nameKey` members,
  helper declaration.
- `src/fastaFromBed/fastaFromBed.cpp` — constructor, `getGffAttribute` helper, `ReportSeq` branch.
- `docs/content/tools/getfasta.rst` — document the flag, the name-only default, and limitations.
- `test/getfasta/test-getfasta.sh` (+ a small sample `.gff`) — see below.

## Testing

Add cases to `test/getfasta/test-getfasta.sh`:

1. Attribute present -> header is the attribute value (`>geneA`).
2. `-nameKey` + `-s` -> strand suffix appended (`>geneA(+)`).
3. Attribute missing on a record -> stderr warning emitted and header falls back to feature type.
4. Non-GFF (BED) input with `-nameKey` -> warning + fallback, no crash.
5. `-nameKey` combined with `-name`/`-nameOnly` -> error exit (mutual exclusivity).

## Out of scope

- GTF `key "value";` attribute syntax.
- URL/percent-decoding of GFF3 attribute values.
- Changing the default GFF behavior of `-name` / `-nameOnly`.
