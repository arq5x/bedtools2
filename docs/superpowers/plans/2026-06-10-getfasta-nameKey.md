# getfasta `-nameKey` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in `getfasta -nameKey <attr>` flag that builds FASTA headers from a named GFF3 column-9 attribute (e.g. `Name`, `ID`), without changing any existing flag's behavior.

**Architecture:** New flag parsed in `fastaFromBedMain.cpp`, plumbed through the `Bed2Fa` constructor as `_useNameKey` / `_nameKey`. A free-standing helper parses the GFF3 attribute string; `ReportSeq` uses it as the name source (name-only output) when the input is GFF, warning + falling back to the feature-type name when the attribute is absent. `-nameKey` is mutually exclusive with `-name`/`-name+`/`-nameOnly` (error + exit).

**Tech Stack:** C++ (g++/clang via top-level `make`), existing shell test harness at `test/getfasta/test-getfasta.sh`.

**Spec:** `docs/superpowers/specs/2026-06-10-getfasta-nameKey-design.md`

---

## File Structure

- `src/fastaFromBed/fastaFromBedMain.cpp` — parse `-nameKey`, mutual-exclusivity check, help text, pass new args to constructor.
- `src/fastaFromBed/fastaFromBed.h` — constructor signature + `_useNameKey` / `_nameKey` members + `getGffAttribute` declaration.
- `src/fastaFromBed/fastaFromBed.cpp` — constructor init, `getGffAttribute` helper, `ReportSeq` branch.
- `src/utils/bedFile/bedFile.h` — add public `isGff()` accessor (currently `_isGff` is private with no getter).
- `test/getfasta/attr.gff` — new sample GFF3 fixture (created in Task 1).
- `test/getfasta/test-getfasta.sh` — new test cases appended before the final `exit` line.
- `docs/content/tools/getfasta.rst` — document the flag.

**Build command (used throughout):** from the repo root, `make -j4 2>&1 | tail -5` (first build is slow; incremental rebuilds after editing only `src/fastaFromBed/*` and `bedFile.h` are fast). The binary lands at `./bin/bedtools`.

**Run the getfasta tests:** `cd test/getfasta && BT=../../bin/bedtools bash test-getfasta.sh`

---

## Task 1: Add `-nameKey` CLI flag, plumbing, and mutual-exclusivity check

This task makes the flag parse and reach the constructor, and enforces mutual exclusivity. Extraction logic comes in Task 2 — here `-nameKey` is wired through but `ReportSeq` does not yet use it.

**Files:**
- Create: `test/getfasta/attr.gff`
- Modify: `src/fastaFromBed/fastaFromBed.h`
- Modify: `src/fastaFromBed/fastaFromBed.cpp:17-35` (constructor signature + initializer list)
- Modify: `src/fastaFromBed/fastaFromBedMain.cpp`
- Test: `test/getfasta/test-getfasta.sh`

- [ ] **Step 1: Create the GFF3 test fixture**

Create `test/getfasta/attr.gff` (TAB-separated; record 1 has `Name`, record 2 does not). Coordinates are 1-based GFF and map onto the existing `t.fa` (chr1, 50 bp):

```
chr1	src	gene	1	10	.	+	.	ID=g1;Name=geneA
chr1	src	mRNA	11	20	.	+	.	ID=g2
```

Verify it is tab-delimited:

Run: `cat -A test/getfasta/attr.gff | head -1`
Expected: fields separated by `^I` (tabs), line ends with `$`.

- [ ] **Step 2: Write the failing test for mutual exclusivity**

Append to `test/getfasta/test-getfasta.sh` immediately BEFORE the final `[[ $FAILURES -eq 0 ]] || exit 1;` line at the end of the file:

```bash
# -nameKey is mutually exclusive with -name/-name+/-nameOnly
echo -e "    getfasta.t19...\c"
$BT getfasta -fi t.fa -bed attr.gff -nameKey Name -name > /dev/null 2> obs
RC=$?
if [ "$RC" == "1" ] && grep -q "cannot be combined" obs; then
    echo ok
else
    FAILURES=$(expr $FAILURES + 1);
    echo fail "rc=$RC"
fi
rm -f obs
```

- [ ] **Step 3: Build and run to verify it fails**

Run: `make -j4 2>&1 | tail -3 && cd test/getfasta && BT=../../bin/bedtools bash test-getfasta.sh 2>&1 | tail -3; cd ../..`
Expected: `getfasta.t19` reports `fail` (today `-nameKey` is an unrecognized parameter, which prints help to stderr and exits 1, but the message is the help text, not "cannot be combined"; the `grep` fails). This confirms the test is meaningful.

- [ ] **Step 4: Add the constructor parameters and members in the header**

In `src/fastaFromBed/fastaFromBed.h`, change the constructor declaration (lines 33-39) to add two parameters at the end (before `isRNA` to keep related name flags together):

```cpp
    Bed2Fa(const string &dbFile, 
           const string &bedFile, const string &fastaOutFile,
           bool useFasta, bool useStrand, 
           bool useBlocks, bool useFullHeader,
           bool useBedOut, bool useName, 
           bool useNamePlus, bool useNameOnly,
           bool useNameKey, const string &nameKey,
           bool isRNA);
```

Add the private members after `_useNameOnly;` (line 61):

```cpp
    bool _useNameOnly;
    bool _useNameKey;
    string _nameKey;
```

Add the helper declaration in the `public:` section after `void ReportSeq(...)` (line 45):

```cpp
    void ReportSeq(const BED &bed, string &seq);
    // Extract the value of GFF3 column-9 attribute `key` from `bed`.
    // Returns true and sets `value` on success; false if not GFF or key absent.
    bool getGffAttribute(const BED &bed, const string &key, string &value);
```

- [ ] **Step 5: Update the constructor definition**

In `src/fastaFromBed/fastaFromBed.cpp`, update the constructor signature (lines 17-23) and initializer list (lines 24-35) to accept and store the new parameters:

```cpp
Bed2Fa::Bed2Fa(const string &dbFile, 
               const string &bedFile, const string &fastaOutFile,
               bool useFasta, bool useStrand, 
               bool useBlocks, bool useFullHeader,
               bool useBedOut, bool useName, 
               bool useNamePlus, bool useNameOnly,
               bool useNameKey, const string &nameKey,
               bool isRNA) :
    _dbFile(dbFile),
    _bedFile(bedFile),
    _fastaOutFile(fastaOutFile),
    _useFasta(useFasta),
    _useStrand(useStrand),
    _useBlocks(useBlocks),
    _useFullHeader(useFullHeader),
    _useBedOut(useBedOut),
    _useName(useName),
    _useNamePlus(useNamePlus),
    _useNameOnly(useNameOnly),
    _useNameKey(useNameKey),
    _nameKey(nameKey),
    _isRNA(isRNA)
{
```

- [ ] **Step 6: Parse `-nameKey` and enforce mutual exclusivity in main**

In `src/fastaFromBed/fastaFromBedMain.cpp`:

(a) Add config variables next to the other name flags (after `bool useNameOnly = false;`, line 45):

```cpp
    bool useNameOnly = false;
    bool useNameKey = false;
    string nameKey;
```

(b) Add the parse branch. Place it right after the `-nameOnly` branch (after line 101):

```cpp
        else if(PARAMETER_CHECK("-nameOnly", 9, parameterLength)) {
            useNameOnly = true;
        }
        else if(PARAMETER_CHECK("-nameKey", 8, parameterLength)) {
            if ((i+1) < argc) {
                useNameKey = true;
                nameKey = argv[i + 1];
                i++;
            }
        }
```

(c) Add the mutual-exclusivity check after the `haveFastaDb`/`haveBed` check (after line 132, before the `if (!haveFastaOut)` block):

```cpp
    if (useNameKey && (useName || useNamePlus || useNameOnly)) {
        cerr << "*****ERROR: -nameKey cannot be combined with -name, "
             << "-name+, or -nameOnly. Choose one. *****" << endl << endl;
        fastafrombed_help();   // prints usage and exit(1)
    }
```

(d) Pass the new args to the constructor (lines 140-146):

```cpp
        Bed2Fa *b2f = new Bed2Fa(fastaDbFile, 
                                 bedFile, fastaOutFile,
                                 useFasta, useStrand, 
                                 useBlocks, useFullHeader,
                                 useBedOut, useName, 
                                 useNamePlus, useNameOnly,
                                 useNameKey, nameKey,
                                 isRNA);
```

(e) Add a help line after the `-nameOnly` line (line 171):

```cpp
    cerr << "\t-nameOnly\tUse the name field for the FASTA header" << endl;
    cerr << "\t-nameKey\tUse the value of the named GFF3 attribute (col 9,"
         << "\n\t\t\te.g. -nameKey Name) for the FASTA header." << endl;
```

- [ ] **Step 7: Add a temporary stub for `getGffAttribute` so it links**

Task 2 fills this in. For now add a minimal stub at the end of `src/fastaFromBed/fastaFromBed.cpp` so the build links:

```cpp
bool Bed2Fa::getGffAttribute(const BED &bed, const string &key, string &value) {
    return false;
}
```

- [ ] **Step 8: Build and run to verify the mutual-exclusivity test passes**

Run: `make -j4 2>&1 | tail -3 && cd test/getfasta && BT=../../bin/bedtools bash test-getfasta.sh 2>&1 | tail -3; cd ../..`
Expected: `getfasta.t19` reports `ok`. All previously-passing tests still `ok`.

- [ ] **Step 9: Commit**

```bash
git add src/fastaFromBed/fastaFromBed.h src/fastaFromBed/fastaFromBed.cpp \
        src/fastaFromBed/fastaFromBedMain.cpp \
        test/getfasta/attr.gff test/getfasta/test-getfasta.sh
git commit -m "getfasta: add -nameKey flag plumbing + mutual-exclusivity check (refs #182)"
```

---

## Task 2: Implement attribute extraction and header logic

**Files:**
- Modify: `src/utils/bedFile/bedFile.h` (add `isGff()` accessor)
- Modify: `src/fastaFromBed/fastaFromBed.cpp` (`getGffAttribute` body + `ReportSeq` branch)
- Test: `test/getfasta/test-getfasta.sh`

- [ ] **Step 1: Write the failing tests for extraction, fallback, and strand**

Append to `test/getfasta/test-getfasta.sh` immediately BEFORE the final `[[ $FAILURES -eq 0 ]] || exit 1;` line (after the t19 block):

```bash
# -nameKey extracts the attribute value (name-only header)
echo -e "    getfasta.t20...\c"
echo ">geneA
aggggggggg" > exp
$BT getfasta -fi t.fa -bed attr.gff -nameKey Name 2>/dev/null | head -2 > obs
check obs exp
rm -f exp obs

# -nameKey with a missing attribute warns and falls back to the feature type
echo -e "    getfasta.t21...\c"
$BT getfasta -fi t.fa -bed attr.gff -nameKey Name 2> err > obs
# record 2 has no Name= ; header falls back to the GFF type "mRNA"
if grep -q "^>mRNA$" obs && grep -q "not found" err; then
    echo ok
else
    FAILURES=$(expr $FAILURES + 1);
    echo fail
fi
rm -f obs err

# -nameKey composes with -s (strand suffix appended)
echo -e "    getfasta.t22...\c"
echo ">geneA(+)" > exp
$BT getfasta -fi t.fa -bed attr.gff -nameKey Name -s 2>/dev/null | head -1 > obs
check obs exp
rm -f exp obs
```

- [ ] **Step 2: Build and run to verify they fail**

Run: `make -j4 2>&1 | tail -3 && cd test/getfasta && BT=../../bin/bedtools bash test-getfasta.sh 2>&1 | tail -5; cd ../..`
Expected: `t20`, `t21`, `t22` report `fail` (the stub returns false, so headers are still feature-type/coords, and no "not found" warning is emitted).

- [ ] **Step 3: Add the public `isGff()` accessor to BedFile**

In `src/utils/bedFile/bedFile.h`, in the `public:` section near `_status`/`_lineNum` (around line 498-499), add:

```cpp
    BedLineStatus _status;
    int _lineNum;
    bool isGff(void) const { return _isGff; }
```

- [ ] **Step 4: Implement `getGffAttribute`**

Replace the stub at the end of `src/fastaFromBed/fastaFromBed.cpp` with the real parser. GFF3 syntax only (`key=value;`), case-sensitive key match, strips one pair of surrounding double-quotes:

```cpp
bool Bed2Fa::getGffAttribute(const BED &bed, const string &key, string &value) {
    // Only GFF records carry a column-9 attribute string.
    if (!_bed->isGff() || bed.fields.size() < 9)
        return false;

    const string &attrs = bed.fields[8];
    size_t pos = 0;
    while (pos < attrs.size()) {
        size_t semi = attrs.find(';', pos);
        string token = attrs.substr(pos, (semi == string::npos)
                                         ? string::npos : semi - pos);
        // trim leading/trailing whitespace
        size_t b = token.find_first_not_of(" \t");
        size_t e = token.find_last_not_of(" \t");
        if (b != string::npos) {
            token = token.substr(b, e - b + 1);
            size_t eq = token.find('=');
            if (eq != string::npos && token.substr(0, eq) == key) {
                string v = token.substr(eq + 1);
                // strip one pair of surrounding double-quotes
                if (v.size() >= 2 && v.front() == '"' && v.back() == '"')
                    v = v.substr(1, v.size() - 2);
                value = v;
                return true;
            }
        }
        if (semi == string::npos) break;
        pos = semi + 1;
    }
    return false;
}
```

- [ ] **Step 5: Add the `_useNameKey` branch to `ReportSeq`**

In `src/fastaFromBed/fastaFromBed.cpp`, replace the header-building `if/else` chain inside `ReportSeq` (lines 80-94, the `else` branch that builds `header`) so the `_useNameKey` case is handled first:

```cpp
        ostringstream header;
        if (_useNameKey)
        {
            string attrVal;
            if (getGffAttribute(bed, _nameKey, attrVal)) {
                header << attrVal;
            }
            else {
                cerr << "WARNING: attribute '" << _nameKey
                     << "' not found for feature at " << bed.chrom << ":"
                     << bed.start << "-" << bed.end
                     << "; using feature name." << endl;
                header << bed.name;
            }
        }
        else if (_useName || _useNamePlus)
        {
            header << bed.name << "::" << bed.chrom << ":" 
                   << bed.start << "-" << bed.end;
        }
        else if (_useNameOnly)
        {
            header << bed.name;
        }
        else 
        {
            header << bed.chrom << ":" 
                   << bed.start << "-" << bed.end;
        }
```

Leave the subsequent `if (_useStrand) { header << "(" << bed.strand << ")"; }` block unchanged — it appends the strand suffix to whatever name source was chosen, so `-nameKey -s` yields `geneA(+)`.

- [ ] **Step 6: Build and run to verify the tests pass**

Run: `make -j4 2>&1 | tail -3 && cd test/getfasta && BT=../../bin/bedtools bash test-getfasta.sh 2>&1 | tail -8; cd ../..`
Expected: `t20`, `t21`, `t22` all `ok`; t19 still `ok`; all earlier tests still `ok`.

- [ ] **Step 7: Commit**

```bash
git add src/utils/bedFile/bedFile.h src/fastaFromBed/fastaFromBed.cpp \
        test/getfasta/test-getfasta.sh
git commit -m "getfasta: implement -nameKey GFF3 attribute extraction (closes #182)"
```

---

## Task 3: Document `-nameKey`

**Files:**
- Modify: `docs/content/tools/getfasta.rst`

- [ ] **Step 1: Find the options/usage section to update**

Run: `grep -n "name\|nameOnly\|Option\|====" docs/content/tools/getfasta.rst | head -30`
Expected: locate the section listing `-name`, `-nameOnly`, etc. (an options table or definition list).

- [ ] **Step 2: Add the `-nameKey` entry**

Add an entry alongside the existing `-name` / `-nameOnly` entries, matching the surrounding markup style. Use this content:

```
**-nameKey**   Use the value of a named GFF3 column-9 attribute as the FASTA
header. Example: ``-nameKey Name`` produces ``>geneA``. The header is the
attribute value only (no coordinates). If the attribute is absent on a
feature, a warning is written to stderr and the header falls back to the
feature-type name. GFF3 ``key=value;`` syntax only (GTF ``key "value";`` is
not parsed). Cannot be combined with ``-name``, ``-name+``, or ``-nameOnly``.
```

- [ ] **Step 3: Commit**

```bash
git add docs/content/tools/getfasta.rst
git commit -m "docs: document getfasta -nameKey (refs #182)"
```

---

## Self-Review Notes

- **Spec coverage:** CLI flag (Task 1 Step 6), name-only default (Task 2 Step 5 `header << attrVal`), GFF3-only parse (Task 2 Step 4), warn+fallback (Task 2 Step 5), mutual exclusivity (Task 1 Step 6c + test t19), `-s` composition (test t22), docs (Task 3), tests (t19-t22). All covered.
- **Type consistency:** `getGffAttribute(const BED&, const string&, string&)` declared in Task 1 Step 4, stubbed in Task 1 Step 7, defined in Task 2 Step 4 — identical signature. `isGff()` added Task 2 Step 3, used Task 2 Step 4. Constructor arg order `useNameKey, nameKey` consistent across header (T1 S4), definition (T1 S5), and call site (T1 S6d).
- **Fixture/sequence check:** `attr.gff` record 1 = chr1:0-10 (after GFF 1-based→0-based), which is `t.fa` line 1 `aggggggggg`; record 2 = chr1:10-20 = line 2, type `mRNA`, no `Name` → fallback `>mRNA`. Matches tests t20/t21.
