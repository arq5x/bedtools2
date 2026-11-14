. ./helpers-intersect.sh formats
set -e

# Input formats, headers, compression, stdin, line endings, buffering and
# multiple databases for bedtools intersect. Expected outputs are written by
# hand or come from the short awk reference implementations below; none is
# produced by bedtools.

N=0
# tc <options> <input class> <derivation>: announce the next test id and, when
# TARGETS_OUT is set, append a row to that coverage table.
tc()
{
	N=$((N + 1))
	ID=$(printf '%03d' $N)
	begin $ID
	if [ -n "${TARGETS_OUT-}" ]; then
		printf '%s\t%s\t%s\t%s\n' "$ID" "$1" "$2" "$3" >> "$TARGETS_OUT"
	fi
}

# tb: turn '|' into TAB so records with spaces in a field stay readable.
tb()
{
	tr '|' '\t'
}

# run <args>: bedtools intersect with stdout in $TMPD/o, stderr in $TMPD/o.err.
run()
{
	RC=0
	"$BT" intersect "$@" > "$TMPD/o" 2> "$TMPD/o.err" || RC=$?
}

# ck: compare $TMPD/o with $TMPD/e; the command must also have exited 0.
ck()
{
	if [ "${RC-0}" -ne 0 ]; then
		FAILURES=$((FAILURES + 1))
		echo "fail (exit $RC)"
		sed 's/^/      /' "$TMPD/o.err"
		RC=0
		return
	fi
	check "$TMPD/o" "$TMPD/e"
}

# cksorted: like ck, but compares $TMPD/o.sorted.
cksorted()
{
	cp "$TMPD/o" "$TMPD/o.raw"
	cp "$TMPD/o.sorted" "$TMPD/o"
	ck
}

# rec <type> <name> <chr> <start0> <end>: one data line of the given type for
# the 0-based half-open interval start0..end. GFF and VCF are written 1-based.
rec()
{
	local t=$1 name=$2 chr=$3 s=$4 e=$5 ref
	case $t in
	bed3) printf '%s\t%s\t%s\n' $chr $s $e ;;
	bed4) printf '%s\t%s\t%s\t%s\n' $chr $s $e $name ;;
	bed5) printf '%s\t%s\t%s\t%s\t10\n' $chr $s $e $name ;;
	bed6) printf '%s\t%s\t%s\t%s\t5\t+\n' $chr $s $e $name ;;
	bed12) printf '%s\t%s\t%s\t%s\t5\t+\t%s\t%s\t0,0,0\t1\t%s\t0\n' $chr $s $e $name $s $e $((e - s)) ;;
	bedplus) printf '%s\t%s\t%s\t%s\tx\ty\tz\n' $chr $s $e $name ;;
	bed6p) printf '%s\t%s\t%s\t%s\t5\t+\tx\ty\n' $chr $s $e $name ;;
	bedgraph) printf '%s\t%s\t%s\t1.5\n' $chr $s $e ;;
	gff) printf '%s\tsrc\tfeat\t%s\t%s\t.\t+\t.\tID=%s;Note=has space\n' $chr $((s + 1)) $e $name ;;
	vcf)
		ref=$(printf '%*s' $((e - s)) '' | tr ' ' N)
		printf '%s\t%s\t%s\t%s\tA\t50\tPASS\t.\n' $chr $((s + 1)) $name $ref ;;
	esac
}

# nul <type>: the NULL record printed for a database of that type.
nul()
{
	case $1 in
	bed3) echo '.|-1|-1' ;;
	bed4) echo '.|-1|-1|.' ;;
	bed5) echo '.|-1|-1|.|-1' ;;
	bed6) echo '.|-1|-1|.|-1|.' ;;
	bed12) echo '.|-1|-1|.|-1|.|.|.|.|.|.|.' ;;
	bedplus) echo '.|-1|-1|.|.|.|.' ;;
	bed6p) echo '.|-1|-1|.|.|.|.|.' ;;
	bedgraph) echo '.|-1|-1|.' ;;
	gff) echo '.|.|.|-1|-1|.|.|.|.' ;;
	vcf) echo '.|-1|.|.|.|.|.|.' ;;
	esac | tb
}

# mkfile <type> <file> <name prefix> <start:end>...: records on chr1.
mkfile()
{
	local t=$1 f=$2 p=$3 i=0 se
	shift 3
	: > "$f"
	if [ "$t" = vcf ]; then
		printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > "$f"
	fi
	for se in "$@"; do
		i=$((i + 1))
		rec $t $p$i chr1 ${se%:*} ${se#*:} >> "$f"
	done
}

# trim <file> <start> <end>: columns 2 and 3 of every line replaced.
trim()
{
	awk -F'\t' -v OFS='\t' -v s=$2 -v e=$3 '{ $2 = s; $3 = e; print }' "$1"
}

# mexp <A type> <B type> <mode>: expected output for the A intervals
# 10-20, 30-40, 50-60 against the B intervals 15-25, 35-38, 70-80. A1 meets
# B1 over 15-20 (5 bases), A2 meets B2 over 35-38 (3 bases), A3 meets nothing.
mexp()
{
	local ta=$1 tb=$2 mode=$3 A1 A2 A3 B1 B2 NUL
	A1=$(rec $ta a1 chr1 10 20)
	A2=$(rec $ta a2 chr1 30 40)
	A3=$(rec $ta a3 chr1 50 60)
	B1=$(rec $tb b1 chr1 15 25)
	B2=$(rec $tb b2 chr1 35 38)
	NUL=$(nul $tb)
	case $mode in
	def)
		{ echo "$A1" | trim /dev/stdin 15 20; echo "$A2" | trim /dev/stdin 35 38; } ;;
	wb)
		{ printf '%s\n' "$A1" | trim /dev/stdin 15 20 | tr '\n' '\t'; echo "$B1"
		  printf '%s\n' "$A2" | trim /dev/stdin 35 38 | tr '\n' '\t'; echo "$B2"; } ;;
	wa) printf '%s\n%s\n' "$A1" "$A2" ;;
	wawb) printf '%s\t%s\n%s\t%s\n' "$A1" "$B1" "$A2" "$B2" ;;
	loj) printf '%s\t%s\n%s\t%s\n%s\t%s\n' "$A1" "$B1" "$A2" "$B2" "$A3" "$NUL" ;;
	wao) printf '%s\t%s\t5\n%s\t%s\t3\n%s\t%s\t0\n' "$A1" "$B1" "$A2" "$B2" "$A3" "$NUL" ;;
	u) printf '%s\n%s\n' "$A1" "$A2" ;;
	v) printf '%s\n' "$A3" ;;
	c) printf '%s\t1\n%s\t1\n%s\t0\n' "$A1" "$A2" "$A3" ;;
	esac > "$TMPD/e"
}

TYPES="bed3 bed4 bed5 bed6 bed12 bedplus bed6p bedgraph gff vcf"
BEDTYPES="bed3 bed4 bed5 bed6 bed12 bedplus bed6p bedgraph"

# ---------------------------------------------------------------------
# Record types: same type in A and B (10-20 etc. as described at mexp).
# ---------------------------------------------------------------------

# Each type as A and B: -wa -wb prints the full A record then the full B record.
for t in $TYPES; do
	tc "-wa -wb" "$t vs $t" "hand"
	mkfile $t "$TMPD/a" a 10:20 30:40 50:60
	mkfile $t "$TMPD/b" b 15:25 35:38 70:80
	mexp $t $t wawb
	run -a "$TMPD/a" -b "$TMPD/b" -wa -wb
	ck
done

# Each type as A and B: -loj appends the NULL record of that type for A3.
for t in $TYPES; do
	tc "-loj" "$t vs $t" "hand"
	mkfile $t "$TMPD/a" a 10:20 30:40 50:60
	mkfile $t "$TMPD/b" b 15:25 35:38 70:80
	mexp $t $t loj
	run -a "$TMPD/a" -b "$TMPD/b" -loj
	ck
done

# Each type as A and B: -wao appends the overlap length (5, 3, 0).
for t in $TYPES; do
	tc "-wao" "$t vs $t" "hand"
	mkfile $t "$TMPD/a" a 10:20 30:40 50:60
	mkfile $t "$TMPD/b" b 15:25 35:38 70:80
	mexp $t $t wao
	run -a "$TMPD/a" -b "$TMPD/b" -wao
	ck
done

# BED-like A types, no -wa/-wb: A is printed with the overlap coordinates and its other columns kept.
for t in $BEDTYPES; do
	tc "(none)" "$t vs $t" "hand: awk replaces columns 2,3 of the A line"
	mkfile $t "$TMPD/a" a 10:20 30:40 50:60
	mkfile $t "$TMPD/b" b 15:25 35:38 70:80
	mexp $t $t def
	run -a "$TMPD/a" -b "$TMPD/b"
	ck
done

# Three BED-like A types with -wb: A with overlap coordinates followed by the B record.
for t in bed3 bed6 bed12; do
	tc "-wb" "$t vs $t" "hand: awk replaces columns 2,3 of the A line"
	mkfile $t "$TMPD/a" a 10:20 30:40 50:60
	mkfile $t "$TMPD/b" b 15:25 35:38 70:80
	mexp $t $t wb
	run -a "$TMPD/a" -b "$TMPD/b" -wb
	ck
done

# -u, -v and -c on a selection of types: A printed once, only A3, or with a count.
for t in bed3 bed6 gff vcf; do
	for m in u v c; do
		tc "-$m" "$t vs $t" "hand"
		mkfile $t "$TMPD/a" a 10:20 30:40 50:60
		mkfile $t "$TMPD/b" b 15:25 35:38 70:80
		mexp $t $t $m
		run -a "$TMPD/a" -b "$TMPD/b" -$m
		ck
	done
done

# Mixed: BED3 query against every other type, -loj prints the NULL of the B type.
for t in bed4 bed5 bed6 bed12 bedplus bed6p bedgraph gff vcf; do
	tc "-loj" "bed3 vs $t" "hand"
	mkfile bed3 "$TMPD/a" a 10:20 30:40 50:60
	mkfile $t "$TMPD/b" b 15:25 35:38 70:80
	mexp bed3 $t loj
	run -a "$TMPD/a" -b "$TMPD/b" -loj
	ck
done

# Mixed: BED3 query against BED12, GFF and VCF, -wao.
for t in bed12 gff vcf; do
	tc "-wao" "bed3 vs $t" "hand"
	mkfile bed3 "$TMPD/a" a 10:20 30:40 50:60
	mkfile $t "$TMPD/b" b 15:25 35:38 70:80
	mexp bed3 $t wao
	run -a "$TMPD/a" -b "$TMPD/b" -wao
	ck
done

# Mixed: every other type as A against a BED3 database, -wa -wb.
for t in bed4 bed5 bed6 bed12 bedplus bed6p bedgraph gff vcf; do
	tc "-wa -wb" "$t vs bed3" "hand"
	mkfile $t "$TMPD/a" a 10:20 30:40 50:60
	mkfile bed3 "$TMPD/b" b 15:25 35:38 70:80
	mexp $t bed3 wawb
	run -a "$TMPD/a" -b "$TMPD/b" -wa -wb
	ck
done

# ---------------------------------------------------------------------
# GFF and GTF: 1-based inclusive coordinates.
# ---------------------------------------------------------------------

tb > "$TMPD/g.gff" <<'EOF'
chr1|src|gene|101|200|.|+|.|ID=g1;Name=my gene
chr1|src|gene|301|400|.|-|.|ID=g2;Name=other gene
EOF
tb > "$TMPD/g.gtf" <<'EOF'
chr1|hav|exon|101|200|.|+|.|gene_id "g1"; transcript_id "t1";
chr1|hav|exon|301|400|.|-|.|gene_id "g2"; transcript_id "t2";
EOF

# GFF A against BED3 B, no flags: start and end of the overlap are printed 1-based (1-based 151..200 is 0-based 150..200).
tc "(none)" "gff vs bed3" "hand"
printf 'chr1\t150\t250\n' > "$TMPD/b"
run -a "$TMPD/g.gff" -b "$TMPD/b"
tb > "$TMPD/e" <<'EOF'
chr1|src|gene|151|200|.|+|.|ID=g1;Name=my gene
EOF
ck

# GTF A against BED3 B, no flags: attribute text with spaces and quotes is kept.
tc "(none)" "gtf vs bed3" "hand"
run -a "$TMPD/g.gtf" -b "$TMPD/b"
tb > "$TMPD/e" <<'EOF'
chr1|hav|exon|151|200|.|+|.|gene_id "g1"; transcript_id "t1";
EOF
ck

# GFF A, -wa -wb: the original GFF line, then the BED3 line.
tc "-wa -wb" "gff vs bed3" "hand"
run -a "$TMPD/g.gff" -b "$TMPD/b" -wa -wb
tb > "$TMPD/e" <<'EOF'
chr1|src|gene|101|200|.|+|.|ID=g1;Name=my gene|chr1|150|250
EOF
ck

# GFF A, -wo: 0-based 150..200 is 50 bases.
tc "-wo" "gff vs bed3" "hand"
run -a "$TMPD/g.gff" -b "$TMPD/b" -wo
tb > "$TMPD/e" <<'EOF'
chr1|src|gene|101|200|.|+|.|ID=g1;Name=my gene|chr1|150|250|50
EOF
ck

# GFF A whose last base (1-based 200 = 0-based 199) is the first base of B 199..300: overlap of 1 base, printed as 200..200.
tc "(none)" "gff vs bed3" "hand"
printf 'chr1\t199\t300\n' > "$TMPD/b"
run -a "$TMPD/g.gff" -b "$TMPD/b"
tb > "$TMPD/e" <<'EOF'
chr1|src|gene|200|200|.|+|.|ID=g1;Name=my gene
EOF
ck

# Same pair with -wo: one base.
tc "-wo" "gff vs bed3" "hand"
run -a "$TMPD/g.gff" -b "$TMPD/b" -wo
tb > "$TMPD/e" <<'EOF'
chr1|src|gene|101|200|.|+|.|ID=g1;Name=my gene|chr1|199|300|1
EOF
ck

# B starting at 200 (0-based) only touches the end of GFF 101..200: no overlap.
tc "(none)" "gff vs bed3" "hand"
printf 'chr1\t200\t300\n' > "$TMPD/b"
run -a "$TMPD/g.gff" -b "$TMPD/b"
: > "$TMPD/e"
ck

# BED6 A against GFF B with -wb: A with overlap coordinates, then the GFF line unchanged (1-based 151..250 is 0-based 150..250).
tc "-wb" "bed6 vs gff" "hand"
printf 'chr1\t100\t200\ta\t5\t+\n' > "$TMPD/a"
tb > "$TMPD/b" <<'EOF'
chr1|src|feat|151|250|.|-|.|ID=f1;Note=has space
EOF
run -a "$TMPD/a" -b "$TMPD/b" -wb
tb > "$TMPD/e" <<'EOF'
chr1|150|200|a|5|+|chr1|src|feat|151|250|.|-|.|ID=f1;Note=has space
EOF
ck

# Same pair with -wo: 50 bases.
tc "-wo" "bed6 vs gff" "hand"
run -a "$TMPD/a" -b "$TMPD/b" -wo
tb > "$TMPD/e" <<'EOF'
chr1|100|200|a|5|+|chr1|src|feat|151|250|.|-|.|ID=f1;Note=has space|50
EOF
ck

# GFF B with a GTF A, both -wa -wb: both lines unchanged.
tc "-wa -wb" "gtf vs gff" "hand"
run -a "$TMPD/g.gtf" -b "$TMPD/g.gff" -wa -wb
tb > "$TMPD/e" <<'EOF'
chr1|hav|exon|101|200|.|+|.|gene_id "g1"; transcript_id "t1";|chr1|src|gene|101|200|.|+|.|ID=g1;Name=my gene
chr1|hav|exon|301|400|.|-|.|gene_id "g2"; transcript_id "t2";|chr1|src|gene|301|400|.|-|.|ID=g2;Name=other gene
EOF
ck

# GFF B with 9 columns, -loj for an A that misses it: NULL GFF record is 9 fields.
tc "-loj" "bed3 vs gff" "hand"
printf 'chr1\t0\t50\n' > "$TMPD/a"
run -a "$TMPD/a" -b "$TMPD/g.gff" -loj
tb > "$TMPD/e" <<'EOF'
chr1|0|50|.|.|.|-1|-1|.|.|.|.
EOF
ck

# GFF B with 8 columns (no attribute column), -loj: NULL GFF record is 8 fields.
tc "-loj" "bed3 vs gff8" "hand"
tb > "$TMPD/b" <<'EOF'
chr1|src|gene|101|200|.|+|.
EOF
run -a "$TMPD/a" -b "$TMPD/b" -loj
tb > "$TMPD/e" <<'EOF'
chr1|0|50|.|.|.|-1|-1|.|.|.
EOF
ck

# GFF A with 8 columns, -wa -wb against BED3.
tc "-wa -wb" "gff8 vs bed3" "hand"
printf 'chr1\t150\t250\n' > "$TMPD/a"
run -a "$TMPD/b" -b "$TMPD/a" -wa -wb
tb > "$TMPD/e" <<'EOF'
chr1|src|gene|101|200|.|+|.|chr1|150|250
EOF
ck

# ---------------------------------------------------------------------
# VCF: end of a variant from REF length, END= or SVLEN=.
# ---------------------------------------------------------------------

tb > "$TMPD/v.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO
chr1|150|snv1|A|G|50|PASS|.
chr1|160|del1|ACGT|A|50|PASS|DP=5
chr1|1000|sv1|N|<DEL>|.|PASS|SVTYPE=DEL;END=1500
chr1|2000|sv2|N|<DEL>|.|PASS|SVTYPE=DEL;SVLEN=-500
chr1|4000|sv3|N|<DEL>|.|PASS|SVTYPE=DEL;SVLEN=-300,-500
EOF
tb > "$TMPD/v.data" <<'EOF'
chr1|150|snv1|A|G|50|PASS|.
chr1|160|del1|ACGT|A|50|PASS|DP=5
chr1|1000|sv1|N|<DEL>|.|PASS|SVTYPE=DEL;END=1500
chr1|2000|sv2|N|<DEL>|.|PASS|SVTYPE=DEL;SVLEN=-500
chr1|4000|sv3|N|<DEL>|.|PASS|SVTYPE=DEL;SVLEN=-300,-500
EOF
# SNV at POS 150 is 0-based 149..150; ACGT at POS 160 is 159..163; END=1500 at POS 1000 is 999..1500.
# The SVLEN variants (sv2, sv3) are only probed far from their end, because bedtools ends them
# at POS-1+|SVLEN|; sv2 spans roughly 1999..2500, sv3 roughly 3999..4500.
# B intervals one base inside the last base of snv1, del1 and sv1:
printf 'chr1\t149\t150\nchr1\t162\t200\nchr1\t1499\t1600\n' > "$TMPD/vhit.bed"
# and starting at the first base after each:
printf 'chr1\t150\t151\nchr1\t163\t200\nchr1\t1500\t1600\n' > "$TMPD/vmiss.bed"
# B well inside sv2, and inside sv3 only if the larger SVLEN (500, not 300) is used:
printf 'chr1\t2100\t2200\nchr1\t4350\t4400\n' > "$TMPD/vsvin.bed"
# B well past both:
printf 'chr1\t2800\t2900\nchr1\t4800\t4900\n' > "$TMPD/vsvout.bed"

# VCF A, -u, B ends inside snv1, del1 and sv1 (END=): those three variants printed as original records, header not printed.
tc "-u" "vcf vs bed3" "hand"
run -a "$TMPD/v.vcf" -b "$TMPD/vhit.bed" -u
head -3 "$TMPD/v.data" > "$TMPD/e"
ck

# VCF A, -u, B starts at the first base after snv1, del1 (REF length) and sv1 (END=): nothing.
tc "-u" "vcf vs bed3" "hand"
run -a "$TMPD/v.vcf" -b "$TMPD/vmiss.bed" -u
: > "$TMPD/e"
ck

# VCF A, -v with the missing B: all five variants (the SVLEN ones are far away).
tc "-v" "vcf vs bed3" "hand"
run -a "$TMPD/v.vcf" -b "$TMPD/vmiss.bed" -v
cp "$TMPD/v.data" "$TMPD/e"
ck

# VCF A, -v with the hitting B: only the two SVLEN variants remain.
tc "-v" "vcf vs bed3" "hand"
run -a "$TMPD/v.vcf" -b "$TMPD/vhit.bed" -v
tail -2 "$TMPD/v.data" > "$TMPD/e"
ck

# SVLEN variants, B inside them: both printed (sv3 only if the larger of SVLEN=-300,-500 is used).
tc "-u" "vcf SVLEN vs bed3" "hand"
run -a "$TMPD/v.vcf" -b "$TMPD/vsvin.bed" -u
tail -2 "$TMPD/v.data" > "$TMPD/e"
ck

# SVLEN variants, B well past both ends: nothing.
tc "-u" "vcf SVLEN vs bed3" "hand"
run -a "$TMPD/v.vcf" -b "$TMPD/vsvout.bed" -u
: > "$TMPD/e"
ck

# VCF A, no flags: the whole original VCF line is printed, not trimmed to the overlap.
tc "(none)" "vcf vs bed3" "hand"
printf 'chr1\t100\t200\n' > "$TMPD/b"
run -a "$TMPD/v.vcf" -b "$TMPD/b"
tb > "$TMPD/e" <<'EOF'
chr1|150|snv1|A|G|50|PASS|.
chr1|160|del1|ACGT|A|50|PASS|DP=5
EOF
ck

# VCF A, -wo: overlap bases are 1 (SNV), 4 (indel by REF length) and 201 (END=1500 variant 999..1500 with B 900..1200).
tc "-wo" "vcf vs bed3" "hand"
printf 'chr1\t100\t200\nchr1\t900\t1200\n' > "$TMPD/b"
run -a "$TMPD/v.vcf" -b "$TMPD/b" -wo
tb > "$TMPD/e" <<'EOF'
chr1|150|snv1|A|G|50|PASS|.|chr1|100|200|1
chr1|160|del1|ACGT|A|50|PASS|DP=5|chr1|100|200|4
chr1|1000|sv1|N|<DEL>|.|PASS|SVTYPE=DEL;END=1500|chr1|900|1200|201
EOF
ck

# VCF A, -header: all three header lines, then the hits (only the first two variants overlap 100..200).
tc "-header" "vcf vs bed3" "hand"
printf 'chr1\t100\t200\n' > "$TMPD/b"
run -a "$TMPD/v.vcf" -b "$TMPD/b" -header
tb > "$TMPD/e" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO
chr1|150|snv1|A|G|50|PASS|.
chr1|160|del1|ACGT|A|50|PASS|DP=5
EOF
ck

# VCF A, -header, B overlaps nothing: the header alone.
tc "-header" "vcf vs bed3" "hand"
printf 'chr9\t1\t2\n' > "$TMPD/b"
run -a "$TMPD/v.vcf" -b "$TMPD/b" -header
head -3 "$TMPD/v.vcf" > "$TMPD/e"
ck

# Header-only VCF as A with -header: the header is still printed.
tc "-header" "header-only vcf vs bed3" "hand"
head -3 "$TMPD/v.vcf" > "$TMPD/ho.vcf"
printf 'chr1\t100\t200\n' > "$TMPD/b"
run -a "$TMPD/ho.vcf" -b "$TMPD/b" -header
cp "$TMPD/ho.vcf" "$TMPD/e"
ck

# Header-only VCF as A without -header: no output.
tc "(none)" "header-only vcf vs bed3" "hand"
run -a "$TMPD/ho.vcf" -b "$TMPD/b"
: > "$TMPD/e"
ck

# BED3 A against a VCF B with -wb: A trimmed to the overlap, then the raw VCF records (header of B never printed).
tc "-wb" "bed3 vs vcf" "hand"
printf 'chr1\t100\t200\n' > "$TMPD/a"
run -a "$TMPD/a" -b "$TMPD/v.vcf" -wb
tb > "$TMPD/e" <<'EOF'
chr1|149|150|chr1|150|snv1|A|G|50|PASS|.
chr1|159|163|chr1|160|del1|ACGT|A|50|PASS|DP=5
EOF
ck

# BED3 A against a VCF B with -wo: 1 base for the SNV, 4 for the 4-base REF.
tc "-wo" "bed3 vs vcf" "hand"
run -a "$TMPD/a" -b "$TMPD/v.vcf" -wo
tb > "$TMPD/e" <<'EOF'
chr1|100|200|chr1|150|snv1|A|G|50|PASS|.|1
chr1|100|200|chr1|160|del1|ACGT|A|50|PASS|DP=5|4
EOF
ck

# VCF A against a VCF B, -wa -wb: a SNV at POS 1100 lies inside END=1500 (999..1500); both printed as in the files.
tc "-wa -wb" "vcf vs vcf" "hand"
tb > "$TMPD/b.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO
chr1|1100|snvb|C|T|10|PASS|.
chr1|3000|snvc|C|T|10|PASS|.
EOF
run -a "$TMPD/v.vcf" -b "$TMPD/b.vcf" -wa -wb
tb > "$TMPD/e" <<'EOF'
chr1|1000|sv1|N|<DEL>|.|PASS|SVTYPE=DEL;END=1500|chr1|1100|snvb|C|T|10|PASS|.
EOF
ck

# VCF with FORMAT and two sample columns as B, -loj: NULL record has 11 fields.
tc "-loj" "bed3 vs vcf with samples" "hand"
tb > "$TMPD/b.vcf" <<'EOF'
##fileformat=VCFv4.2
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO|FORMAT|S1|S2
chr1|1100|snvb|C|T|10|PASS|.|GT|0/1|1/1
EOF
printf 'chr1\t0\t5\n' > "$TMPD/a"
run -a "$TMPD/a" -b "$TMPD/b.vcf" -loj
tb > "$TMPD/e" <<'EOF'
chr1|0|5|.|-1|.|.|.|.|.|.|.|.|.
EOF
ck

# VCF with samples as A against BED3, -wa: the full line with genotype columns.
tc "-wa" "vcf with samples vs bed3" "hand"
printf 'chr1\t1099\t1100\n' > "$TMPD/a"
run -a "$TMPD/b.vcf" -b "$TMPD/a" -wa
tb > "$TMPD/e" <<'EOF'
chr1|1100|snvb|C|T|10|PASS|.|GT|0/1|1/1
EOF
ck

# ---------------------------------------------------------------------
# Headers
# ---------------------------------------------------------------------

tb > "$TMPD/h.bed" <<'EOF'
#comment line
track name=mytrack description="a track"
browser position chr1:1-100
chr1|10|20
#commented out record
chr1|30|40
EOF
printf 'chr1\t15\t35\n' > "$TMPD/b"

# '#', 'track' and 'browser' lines at the top of A are printed with -header.
tc "-header" "bed3 A with headers" "hand"
run -a "$TMPD/h.bed" -b "$TMPD/b" -header
tb > "$TMPD/e" <<'EOF'
#comment line
track name=mytrack description="a track"
browser position chr1:1-100
chr1|15|20
chr1|30|35
EOF
ck

# Without -header none of the header lines nor the mid-file '#' line is printed.
tc "(none)" "bed3 A with headers" "hand"
run -a "$TMPD/h.bed" -b "$TMPD/b"
tb > "$TMPD/e" <<'EOF'
chr1|15|20
chr1|30|35
EOF
ck

# A '#' line after data is a commented-out record, not part of the header: absent even with -header (see the first header test).
tc "-header -wa -wb" "bed3 A with headers" "hand"
run -a "$TMPD/h.bed" -b "$TMPD/b" -header -wa -wb
tb > "$TMPD/e" <<'EOF'
#comment line
track name=mytrack description="a track"
browser position chr1:1-100
chr1|10|20|chr1|15|35
chr1|30|40|chr1|15|35
EOF
ck

# Headers in B are never printed, with or without -header.
tc "-header" "bed3 B with headers" "hand"
run -a "$TMPD/b" -b "$TMPD/h.bed" -header
tb > "$TMPD/e" <<'EOF'
chr1|15|20
chr1|30|35
EOF
ck

# A header in B is skipped when B is read: -wb output has only data columns.
tc "-wb" "bed3 B with headers" "hand"
run -a "$TMPD/b" -b "$TMPD/h.bed" -wb
tb > "$TMPD/e" <<'EOF'
chr1|15|20|chr1|10|20
chr1|30|35|chr1|30|40
EOF
ck

# Headers in both files, -header: only A's.
tc "-header" "headers in A and B" "hand"
printf '#only A\nchr1\t10\t20\n' > "$TMPD/a"
printf '#only B\nchr1\t15\t25\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -header
tb > "$TMPD/e" <<'EOF'
#only A
chr1|15|20
EOF
ck

# -header when nothing overlaps: the header is still printed.
tc "-header" "bed3 A with header, no hits" "hand"
printf 'chr1\t500\t600\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -header
printf '#only A\n' > "$TMPD/e"
ck

# -header with -c: header, then counts.
tc "-header -c" "bed3 A with header" "hand"
printf 'chr1\t15\t25\n' > "$TMPD/b"
printf '#only A\nchr1\t10\t20\nchr1\t100\t110\n' > "$TMPD/a"
run -a "$TMPD/a" -b "$TMPD/b" -header -c
tb > "$TMPD/e" <<'EOF'
#only A
chr1|10|20|1
chr1|100|110|0
EOF
ck

# -header with -loj.
tc "-header -loj" "bed3 A with header" "hand"
run -a "$TMPD/a" -b "$TMPD/b" -header -loj
tb > "$TMPD/e" <<'EOF'
#only A
chr1|10|20|chr1|15|25
chr1|100|110|.|-1|-1
EOF
ck

# -header with -v: header, then A records without overlap.
tc "-header -v" "bed3 A with header" "hand"
run -a "$TMPD/a" -b "$TMPD/b" -header -v
tb > "$TMPD/e" <<'EOF'
#only A
chr1|100|110
EOF
ck

# A file holding only header lines, -header: the header is printed.
tc "-header" "header-only bed A" "hand"
printf '#c1\ntrack name=x\n' > "$TMPD/ho.bed"
run -a "$TMPD/ho.bed" -b "$TMPD/b" -header
printf '#c1\ntrack name=x\n' > "$TMPD/e"
ck

# A file holding only header lines, no -header: nothing.
tc "(none)" "header-only bed A" "hand"
run -a "$TMPD/ho.bed" -b "$TMPD/b"
: > "$TMPD/e"
ck

# Header-only B: no hits, so -u prints nothing.
tc "-u" "header-only bed B" "hand"
run -a "$TMPD/a" -b "$TMPD/ho.bed" -u
: > "$TMPD/e"
ck

# A column-title line starting with "chrom" and a TAB is a header.
tc "-header" "bed3 A with chrom header" "hand"
printf 'chrom\tstart\tend\nchr1\t10\t20\n' > "$TMPD/a"
printf 'chr1\t15\t25\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -header
tb > "$TMPD/e" <<'EOF'
chrom|start|end
chr1|15|20
EOF
ck

# Case-insensitive 'Track' and 'BROWSER' header lines.
tc "-header" "bed3 A with mixed-case headers" "hand"
printf 'Track name=x\nBROWSER hide all\nchr1\t10\t20\n' > "$TMPD/a"
run -a "$TMPD/a" -b "$TMPD/b" -header
printf 'Track name=x\nBROWSER hide all\nchr1\t15\t20\n' > "$TMPD/e"
ck

# A chromosome named like 'chr1' is not mistaken for a 'chr' header line.
tc "(none)" "bed3 chrom names starting with chr" "hand"
printf 'chr1\t10\t20\nchrUn_KI270302v1\t10\t20\n' > "$TMPD/a"
printf 'chr1\t15\t25\nchrUn_KI270302v1\t15\t25\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b"
tb > "$TMPD/e" <<'EOF'
chr1|15|20
chrUn_KI270302v1|15|20
EOF
ck

# ---------------------------------------------------------------------
# Compressed input
# ---------------------------------------------------------------------

printf '#hdr\nchr1\t10\t20\tn1\nchr1\t30\t40\tn2\nchr1\t50\t60\tn3\n' > "$TMPD/ga.bed"
printf 'chr1\t15\t35\tx\n' > "$TMPD/gb.bed"
gzip -c "$TMPD/ga.bed" > "$TMPD/ga.bed.gz"
gzip -c "$TMPD/gb.bed" > "$TMPD/gb.bed.gz"
$htsutil bgzfcompress "$TMPD/ga.bed" "$TMPD/ga.bed.bgz"
$htsutil bgzfcompress "$TMPD/gb.bed" "$TMPD/gb.bed.bgz"
tb > "$TMPD/e.wawb" <<'EOF'
chr1|10|20|n1|chr1|15|35|x
chr1|30|40|n2|chr1|15|35|x
EOF

# gzipped A, plain B.
tc "-wa -wb" "gzip A" "hand"
run -a "$TMPD/ga.bed.gz" -b "$TMPD/gb.bed" -wa -wb
cp "$TMPD/e.wawb" "$TMPD/e"
ck

# plain A, gzipped B.
tc "-wa -wb" "gzip B" "hand"
run -a "$TMPD/ga.bed" -b "$TMPD/gb.bed.gz" -wa -wb
cp "$TMPD/e.wawb" "$TMPD/e"
ck

# both gzipped.
tc "-wa -wb" "gzip A and B" "hand"
run -a "$TMPD/ga.bed.gz" -b "$TMPD/gb.bed.gz" -wa -wb
cp "$TMPD/e.wawb" "$TMPD/e"
ck

# bgzipped A, plain B.
tc "-wa -wb" "bgzip A" "hand"
run -a "$TMPD/ga.bed.bgz" -b "$TMPD/gb.bed" -wa -wb
cp "$TMPD/e.wawb" "$TMPD/e"
ck

# plain A, bgzipped B.
tc "-wa -wb" "bgzip B" "hand"
run -a "$TMPD/ga.bed" -b "$TMPD/gb.bed.bgz" -wa -wb
cp "$TMPD/e.wawb" "$TMPD/e"
ck

# bgzipped A and B.
tc "-wa -wb" "bgzip A and B" "hand"
run -a "$TMPD/ga.bed.bgz" -b "$TMPD/gb.bed.bgz" -wa -wb
cp "$TMPD/e.wawb" "$TMPD/e"
ck

# Two concatenated gzip members as A: records from both members are read.
tc "-wa -wb" "multi-member gzip A" "hand"
printf 'chr1\t10\t20\tn1\n' | gzip -c > "$TMPD/mm.gz"
printf 'chr1\t30\t40\tn2\nchr1\t50\t60\tn3\n' | gzip -c >> "$TMPD/mm.gz"
printf 'chr1\t15\t55\tx\n' > "$TMPD/mmb.bed"
run -a "$TMPD/mm.gz" -b "$TMPD/mmb.bed" -wa -wb
tb > "$TMPD/e" <<'EOT'
chr1|10|20|n1|chr1|15|55|x
chr1|30|40|n2|chr1|15|55|x
chr1|50|60|n3|chr1|15|55|x
EOT
ck

# Two concatenated gzip members as B: hits come from both members.
tc "-wa -wb" "multi-member gzip B" "hand"
printf 'chr1\t0\t100\n' > "$TMPD/mma.bed"
printf 'chr1\t12\t14\tp\n' | gzip -c > "$TMPD/mmb.gz"
printf 'chr1\t30\t32\tq\n' | gzip -c >> "$TMPD/mmb.gz"
run -a "$TMPD/mma.bed" -b "$TMPD/mmb.gz" -wa -wb
tb > "$TMPD/e" <<'EOT'
chr1|0|100|chr1|12|14|p
chr1|0|100|chr1|30|32|q
EOT
ck

# gzipped A with a '#' header and -header: the header is printed.
tc "-header" "gzip A with header" "hand"
run -a "$TMPD/ga.bed.gz" -b "$TMPD/gb.bed" -header
tb > "$TMPD/e" <<'EOF'
#hdr
chr1|15|20|n1
chr1|30|35|n2
EOF
ck

# bgzipped A with a '#' header and -header.
tc "-header" "bgzip A with header" "hand"
run -a "$TMPD/ga.bed.bgz" -b "$TMPD/gb.bed" -header
tb > "$TMPD/e" <<'EOF'
#hdr
chr1|15|20|n1
chr1|30|35|n2
EOF
ck

# gzipped A with a header, no -header: header absent.
tc "(none)" "gzip A with header" "hand"
run -a "$TMPD/ga.bed.gz" -b "$TMPD/gb.bed"
tb > "$TMPD/e" <<'EOF'
chr1|15|20|n1
chr1|30|35|n2
EOF
ck

# gzipped header-only A with -header.
tc "-header" "gzip header-only A" "hand"
printf '#c1\ntrack name=x\n' | gzip -c > "$TMPD/ho.bed.gz"
run -a "$TMPD/ho.bed.gz" -b "$TMPD/gb.bed" -header
printf '#c1\ntrack name=x\n' > "$TMPD/e"
ck

# gzipped VCF A with -header and -u.
tc "-header -u" "gzip vcf A" "hand"
gzip -c "$TMPD/v.vcf" > "$TMPD/v.vcf.gz"
printf 'chr1\t149\t150\n' > "$TMPD/b"
run -a "$TMPD/v.vcf.gz" -b "$TMPD/b" -header -u
tb > "$TMPD/e" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO
chr1|150|snv1|A|G|50|PASS|.
EOF
ck

# gzipped GFF A.
tc "-wa" "gzip gff A" "hand"
gzip -c "$TMPD/g.gtf" > "$TMPD/g.gtf.gz"
printf 'chr1\t150\t250\n' > "$TMPD/b"
run -a "$TMPD/g.gtf.gz" -b "$TMPD/b" -wa
tb > "$TMPD/e" <<'EOF'
chr1|hav|exon|101|200|.|+|.|gene_id "g1"; transcript_id "t1";
EOF
ck

# gzip and bgzip databases together, with a plain one, in order: ids 1, 2, 3.
tc "-wa -wb multiple -b" "plain, gzip, bgzip databases" "hand"
printf 'chr1\t12\t14\td1\n' > "$TMPD/d1.bed"
printf 'chr1\t13\t15\td2\n' | gzip -c > "$TMPD/d2.bed.gz"
printf 'chr1\t14\t16\td3\n' > "$TMPD/d3.src"
$htsutil bgzfcompress "$TMPD/d3.src" "$TMPD/d3.bed.bgz"
printf 'chr1\t10\t20\tn1\n' > "$TMPD/a"
run -a "$TMPD/a" -b "$TMPD/d1.bed" "$TMPD/d2.bed.gz" "$TMPD/d3.bed.bgz" -wa -wb
tb > "$TMPD/e" <<'EOF'
chr1|10|20|n1|1|chr1|12|14|d1
chr1|10|20|n1|2|chr1|13|15|d2
chr1|10|20|n1|3|chr1|14|16|d3
EOF
ck

# ---------------------------------------------------------------------
# stdin
# ---------------------------------------------------------------------

printf 'chr1\t10\t20\nchr1\t30\t40\n' > "$TMPD/a"
printf 'chr1\t15\t35\n' > "$TMPD/b"
tb > "$TMPD/e.in" <<'EOF'
chr1|15|20
chr1|30|35
EOF

# -a stdin.
tc "-a stdin" "stdin A" "hand"
RC=0; "$BT" intersect -a stdin -b "$TMPD/b" < "$TMPD/a" > "$TMPD/o" 2>/dev/null || RC=$?
cp "$TMPD/e.in" "$TMPD/e"
ck

# -a - .
tc "-a -" "stdin A" "hand"
RC=0; "$BT" intersect -a - -b "$TMPD/b" < "$TMPD/a" > "$TMPD/o" 2>/dev/null || RC=$?
cp "$TMPD/e.in" "$TMPD/e"
ck

# -b stdin.
tc "-b stdin" "stdin B" "hand"
RC=0; "$BT" intersect -a "$TMPD/a" -b stdin < "$TMPD/b" > "$TMPD/o" 2>/dev/null || RC=$?
cp "$TMPD/e.in" "$TMPD/e"
ck

# -b - .
tc "-b -" "stdin B" "hand"
RC=0; "$BT" intersect -a "$TMPD/a" -b - < "$TMPD/b" > "$TMPD/o" 2>/dev/null || RC=$?
cp "$TMPD/e.in" "$TMPD/e"
ck

# stdin from a pipe, not a file redirect.
tc "-a stdin" "stdin A from a pipe" "hand"
RC=0; cat "$TMPD/a" | "$BT" intersect -a stdin -b "$TMPD/b" > "$TMPD/o" 2>/dev/null || RC=$?
cp "$TMPD/e.in" "$TMPD/e"
ck

# gzipped data on stdin.
tc "-a stdin" "gzip on stdin" "hand"
RC=0; gzip -c "$TMPD/a" | "$BT" intersect -a stdin -b "$TMPD/b" > "$TMPD/o" 2>/dev/null || RC=$?
cp "$TMPD/e.in" "$TMPD/e"
ck

# stdin with a header and -header.
tc "-a stdin -header" "stdin A with header" "hand"
RC=0; printf '#h\nchr1\t10\t20\n' | "$BT" intersect -a stdin -b "$TMPD/b" -header > "$TMPD/o" 2>/dev/null || RC=$?
tb > "$TMPD/e" <<'EOF'
#h
chr1|15|20
EOF
ck

# VCF on stdin as A, -u: the variant line.
tc "-a stdin -u" "stdin vcf A" "hand"
printf 'chr1\t149\t150\n' > "$TMPD/b"
RC=0; "$BT" intersect -a stdin -b "$TMPD/b" -u < "$TMPD/v.vcf" > "$TMPD/o" 2>/dev/null || RC=$?
printf 'chr1\t150\tsnv1\tA\tG\t50\tPASS\t.\n' > "$TMPD/e"
ck

# stdin as the only database among several, in the middle position of the id list.
tc "-wa -wb -b stdin" "stdin B with other dbs" "hand"
printf 'chr1\t12\t14\n' > "$TMPD/d1.bed"
printf 'chr1\t13\t15\n' > "$TMPD/d2.bed"
printf 'chr1\t10\t20\n' > "$TMPD/a"
RC=0; "$BT" intersect -a "$TMPD/a" -b "$TMPD/d1.bed" stdin -wa -wb < "$TMPD/d2.bed" > "$TMPD/o" 2>/dev/null || RC=$?
tb > "$TMPD/e" <<'EOF'
chr1|10|20|1|chr1|12|14
chr1|10|20|2|chr1|13|15
EOF
ck

# ---------------------------------------------------------------------
# Line endings
# ---------------------------------------------------------------------

# CRLF BED3 in A and B: the CR is not part of the end coordinate and not printed.
tc "(none)" "CRLF bed3" "hand"
printf 'chr1\t10\t20\r\nchr1\t30\t40\r\n' > "$TMPD/a"
printf 'chr1\t15\t35\r\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b"
printf 'chr1\t15\t20\nchr1\t30\t35\n' > "$TMPD/e"
ck

# CRLF BED6: the strand column loses its CR in A and B.
tc "-wa -wb" "CRLF bed6" "hand"
printf 'chr1\t10\t20\tn1\t5\t+\r\n' > "$TMPD/a"
printf 'chr1\t15\t35\tx\t5\t-\r\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wa -wb
printf 'chr1\t10\t20\tn1\t5\t+\tchr1\t15\t35\tx\t5\t-\n' > "$TMPD/e"
ck

# CRLF BED4, -loj, -wao: NULL records and overlap column after a stripped CR.
tc "-wao" "CRLF bed4" "hand"
printf 'chr1\t10\t20\tn1\r\nchr1\t50\t60\tn2\r\n' > "$TMPD/a"
printf 'chr1\t15\t35\tx\r\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wao
tb > "$TMPD/e" <<'EOF'
chr1|10|20|n1|chr1|15|35|x|5
chr1|50|60|n2|.|-1|-1|.|0
EOF
ck

# CRLF GFF A with an attribute column.
tc "-wa" "CRLF gff" "hand"
printf 'chr1\tsrc\tg\t11\t20\t.\t+\t.\tID=a b\r\n' > "$TMPD/a"
printf 'chr1\t15\t35\r\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wa
printf 'chr1\tsrc\tg\t11\t20\t.\t+\t.\tID=a b\n' > "$TMPD/e"
ck

# CRLF VCF with a CRLF header, -header.
tc "-header -u" "CRLF vcf" "hand"
printf '##fileformat=VCFv4.2\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\nchr1\t11\tr\tA\tG\t5\tPASS\t.\r\n' > "$TMPD/a"
printf 'chr1\t10\t11\r\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -header -u
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t11\tr\tA\tG\t5\tPASS\t.\n' > "$TMPD/e"
ck

# CRLF header in A with -header.
tc "-header" "CRLF bed3 with header" "hand"
printf '#h\r\nchr1\t10\t20\r\n' > "$TMPD/a"
printf 'chr1\t15\t35\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -header
printf '#h\nchr1\t15\t20\n' > "$TMPD/e"
ck

# LF in A and CRLF in B, mixed.
tc "-wa -wb" "LF A, CRLF B" "hand"
printf 'chr1\t10\t20\n' > "$TMPD/a"
printf 'chr1\t15\t35\r\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wa -wb
printf 'chr1\t10\t20\tchr1\t15\t35\n' > "$TMPD/e"
ck

# No newline after the last A line.
tc "-wa -wb" "no final newline in A" "hand"
printf 'chr1\t10\t20\nchr1\t30\t40' > "$TMPD/a"
printf 'chr1\t15\t35\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wa -wb
printf 'chr1\t10\t20\tchr1\t15\t35\nchr1\t30\t40\tchr1\t15\t35\n' > "$TMPD/e"
ck

# No newline after the last B line.
tc "-wa -wb" "no final newline in B" "hand"
printf 'chr1\t10\t20\n' > "$TMPD/a"
printf 'chr1\t15\t25\nchr1\t18\t19' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wa -wb
printf 'chr1\t10\t20\tchr1\t15\t25\nchr1\t10\t20\tchr1\t18\t19\n' > "$TMPD/e"
ck

# No newline after the last line of a BED6 A, -c.
tc "-c" "no final newline in bed6" "hand"
printf 'chr1\t10\t20\tn1\t5\t+' > "$TMPD/a"
printf 'chr1\t15\t25\tx\t5\t+' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -c
printf 'chr1\t10\t20\tn1\t5\t+\t1\n' > "$TMPD/e"
ck

# No newline after the last VCF line, A and B.
tc "-wa -wb" "no final newline in vcf" "hand"
printf '##fileformat=VCFv4.2\nchr1\t11\tr\tA\tG\t5\tPASS\t.' > "$TMPD/a"
printf 'chr1\t10\t11' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wa -wb
printf 'chr1\t11\tr\tA\tG\t5\tPASS\t.\tchr1\t10\t11\n' > "$TMPD/e"
ck

# No newline after the last of two GFF lines, gzipped.
tc "-wa" "no final newline in gzip gff" "hand"
printf 'chr1\tsrc\tg\t1\t5\t.\t+\t.\tID=z\nchr1\tsrc\tg\t11\t20\t.\t+\t.\tID=a' | gzip -c > "$TMPD/a.gz"
printf 'chr1\t15\t35' > "$TMPD/b"
run -a "$TMPD/a.gz" -b "$TMPD/b" -wa
printf 'chr1\tsrc\tg\t11\t20\t.\t+\t.\tID=a\n' > "$TMPD/e"
ck

# No newline after a header-only file with -header: header printed with a newline.
tc "-header" "header-only without final newline" "hand"
printf '#c1' > "$TMPD/a"
run -a "$TMPD/a" -b "$TMPD/b" -header
printf '#c1\n' > "$TMPD/e"
ck

# ---------------------------------------------------------------------
# Chromosome names and coordinates
# ---------------------------------------------------------------------

# Names with underscore, hyphen, asterisk and colon, unsorted mode.
tc "-wa -wb" "unusual chromosome names" "hand"
tb > "$TMPD/a" <<'EOF'
chrUn_KI270302v1|10|20
HLA-A*01:01|5|15
chrM|1|10
EOF
tb > "$TMPD/b" <<'EOF'
HLA-A*01:01|10|20
chrUn_KI270302v1|15|25
chrM|20|30
EOF
run -nonamecheck -a "$TMPD/a" -b "$TMPD/b" -wa -wb
tb > "$TMPD/e" <<'EOF'
chrUn_KI270302v1|10|20|chrUn_KI270302v1|15|25
HLA-A*01:01|5|15|HLA-A*01:01|10|20
EOF
ck

# The same names in lexicographic order with -sorted.
tc "-sorted -wa -wb" "unusual chromosome names" "hand"
tb > "$TMPD/a" <<'EOF'
HLA-A*01:01|5|15
chrM|1|10
chrUn_KI270302v1|10|20
EOF
tb > "$TMPD/b" <<'EOF'
HLA-A*01:01|10|20
chrM|20|30
chrUn_KI270302v1|15|25
EOF
run -nonamecheck -sorted -a "$TMPD/a" -b "$TMPD/b" -wa -wb
tb > "$TMPD/e" <<'EOF'
HLA-A*01:01|5|15|HLA-A*01:01|10|20
chrUn_KI270302v1|10|20|chrUn_KI270302v1|15|25
EOF
ck

# -loj on those names prints the NULL record for chrM.
tc "-loj" "unusual chromosome names" "hand"
run -nonamecheck -a "$TMPD/a" -b "$TMPD/b" -loj
tb > "$TMPD/e" <<'EOF'
HLA-A*01:01|5|15|HLA-A*01:01|10|20
chrM|1|10|.|-1|-1
chrUn_KI270302v1|10|20|chrUn_KI270302v1|15|25
EOF
ck

# Purely numeric names match each other; '1' is not 'chr1'.
tc "-wa -wb" "numeric chromosome names" "hand"
tb > "$TMPD/a" <<'EOF'
1|10|20
10|10|20
2|10|20
EOF
tb > "$TMPD/b" <<'EOF'
1|15|25
2|15|25
chr1|15|25
EOF
run -nonamecheck -a "$TMPD/a" -b "$TMPD/b" -wa -wb
tb > "$TMPD/e" <<'EOF'
1|10|20|1|15|25
2|10|20|2|15|25
EOF
ck

# Numeric names in lexicographic order (1, 10, 2) with -sorted.
tc "-sorted -wa -wb" "numeric chromosome names" "hand"
tb > "$TMPD/b" <<'EOF'
1|15|25
10|15|25
2|15|25
EOF
run -nonamecheck -sorted -a "$TMPD/a" -b "$TMPD/b" -wa -wb
tb > "$TMPD/e" <<'EOF'
1|10|20|1|15|25
10|10|20|10|15|25
2|10|20|2|15|25
EOF
ck

# A 'chr1' query against '1' database finds nothing.
tc "-u" "chr1 vs 1 chromosome names" "hand"
printf 'chr1\t10\t20\n' > "$TMPD/a"
printf '1\t10\t20\n' > "$TMPD/b"
run -nonamecheck -a "$TMPD/a" -b "$TMPD/b" -u
: > "$TMPD/e"
ck

# ... and -v reports it.
tc "-v" "chr1 vs 1 chromosome names" "hand"
run -nonamecheck -a "$TMPD/a" -b "$TMPD/b" -v
printf 'chr1\t10\t20\n' > "$TMPD/e"
ck

# Coordinates above 2^31 and 2^32 in unsorted mode.
tc "-wo" "coordinates above 2^31" "hand"
tb > "$TMPD/a" <<'EOF'
chr1|3000000000|3000000100
chr1|4294967290|4294967400
chr1|5000000000|5000000100
EOF
tb > "$TMPD/b" <<'EOF'
chr1|3000000050|3000000150
chr1|4294967300|4294967500
EOF
run -a "$TMPD/a" -b "$TMPD/b" -wo
tb > "$TMPD/e" <<'EOF'
chr1|3000000000|3000000100|chr1|3000000050|3000000150|50
chr1|4294967290|4294967400|chr1|4294967300|4294967500|100
EOF
ck

# The same with -sorted.
tc "-sorted -wo" "coordinates above 2^31" "hand"
run -sorted -a "$TMPD/a" -b "$TMPD/b" -wo
ck

# Default mode prints overlap coordinates above 2^32 in full.
tc "(none)" "coordinates above 2^31" "hand"
run -a "$TMPD/a" -b "$TMPD/b"
tb > "$TMPD/e" <<'EOF'
chr1|3000000050|3000000100
chr1|4294967300|4294967400
EOF
ck

# -c, -v with large coordinates.
tc "-c" "coordinates above 2^31" "hand"
run -a "$TMPD/a" -b "$TMPD/b" -c
tb > "$TMPD/e" <<'EOF'
chr1|3000000000|3000000100|1
chr1|4294967290|4294967400|1
chr1|5000000000|5000000100|0
EOF
ck

# Coordinates near 10^12 in -sorted mode (the unsorted tree cannot hold them).
tc "-sorted -wo" "coordinates near 10^12" "hand"
printf 'chr1\t1000000000000\t1000000000100\n' > "$TMPD/a"
printf 'chr1\t1000000000050\t1000000000150\n' > "$TMPD/b"
run -sorted -a "$TMPD/a" -b "$TMPD/b" -wo
printf 'chr1\t1000000000000\t1000000000100\tchr1\t1000000000050\t1000000000150\t50\n' > "$TMPD/e"
ck

# A GFF record with 1-based start above 2^32.
tc "(none)" "gff coordinates above 2^32" "hand"
printf 'chr1\tsrc\tg\t5000000001\t5000000100\t.\t+\t.\tID=a\n' > "$TMPD/a"
printf 'chr1\t5000000050\t5000000200\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b"
printf 'chr1\tsrc\tg\t5000000051\t5000000100\t.\t+\t.\tID=a\n' > "$TMPD/e"
ck

# ---------------------------------------------------------------------
# Long lines and buffering
# ---------------------------------------------------------------------

LONG=$(awk 'BEGIN { s = ""; for (i = 0; i < 70000; i++) s = s "x"; print s }')

# A 70000-character column in A and B is passed through whole.
tc "-wa -wb" "bed+ with 70000-character column" "hand"
printf 'chr1\t10\t20\tname\t%s\n' "$LONG" > "$TMPD/a"
printf 'chr1\t15\t25\tn2\t%s\n' "$LONG" > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wa -wb
printf 'chr1\t10\t20\tname\t%s\tchr1\t15\t25\tn2\t%s\n' "$LONG" "$LONG" > "$TMPD/e"
ck

# A long line preceded by short lines.
tc "-wa" "bed+ with a later 70000-character line" "hand"
{ printf 'chr1\t1\t2\ts\tx\n'; printf 'chr1\t10\t20\tname\t%s\n' "$LONG"; printf 'chr1\t100\t200\tt\ty\n'; } > "$TMPD/a"
printf 'chr1\t15\t25\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -wa
printf 'chr1\t10\t20\tname\t%s\n' "$LONG" > "$TMPD/e"
ck

# Long line in a BED12 block-size column.
tc "-wa -wb" "bed12 with long blockSizes" "hand"
BS=$(awk 'BEGIN { s = ""; for (i = 0; i < 6000; i++) s = s "1"; print s }')
printf 'chr1\t10\t20\tn\t5\t+\t10\t20\t0,0,0\t1\t%s\t0\n' "$BS" > "$TMPD/a"
run -a "$TMPD/a" -b "$TMPD/a" -wa -wb
printf 'chr1\t10\t20\tn\t5\t+\t10\t20\t0,0,0\t1\t%s\t0\tchr1\t10\t20\tn\t5\t+\t10\t20\t0,0,0\t1\t%s\t0\n' "$BS" "$BS" > "$TMPD/e"
ck

# Long line, -c.
tc "-c" "bed+ with 70000-character column" "hand"
printf 'chr1\t10\t20\tname\t%s\n' "$LONG" > "$TMPD/a"
printf 'chr1\t15\t25\n' > "$TMPD/b"
run -a "$TMPD/a" -b "$TMPD/b" -c
printf 'chr1\t10\t20\tname\t%s\t1\n' "$LONG" > "$TMPD/e"
ck

# Input larger than every buffer size: 60000 A records and 60000 B records, each A meeting
# exactly one B over 25 bases. Expected output is generated with awk from the formula
# A_i = [100i, 100i+50), B_i = [100i+25, 100i+75).
awk 'BEGIN { for (i = 0; i < 60000; i++) printf "chr1\t%d\t%d\n", i * 100, i * 100 + 50 }' > "$TMPD/biga"
awk 'BEGIN { for (i = 0; i < 60000; i++) printf "chr1\t%d\t%d\n", i * 100 + 25, i * 100 + 75 }' > "$TMPD/bigb"
awk 'BEGIN { for (i = 0; i < 60000; i++) printf "chr1\t%d\t%d\tchr1\t%d\t%d\n", i * 100, i * 100 + 50, i * 100 + 25, i * 100 + 75 }' > "$TMPD/bige"

# Buffered output on the large input.
tc "-wa -wb" "60000 records" "reference awk"
run -a "$TMPD/biga" -b "$TMPD/bigb" -wa -wb
cp "$TMPD/bige" "$TMPD/e"
ck

# -nobuf gives the same bytes as the formula.
tc "-nobuf -wa -wb" "60000 records" "reference awk"
run -nobuf -a "$TMPD/biga" -b "$TMPD/bigb" -wa -wb
cp "$TMPD/bige" "$TMPD/e"
ck

# -iobuf values below the 4096-byte scan buffer are raised to it, so 8 behaves like 4096; 4096 and
# above change the plain-text read buffer, and the 1.2 MB input crosses every size listed here.
for sz in 8 4096 64K 1M 1048577; do
	# -iobuf of a given size on the large input.
	tc "-iobuf $sz -wa -wb" "60000 records" "reference awk"
	run -iobuf $sz -a "$TMPD/biga" -b "$TMPD/bigb" -wa -wb
	cp "$TMPD/bige" "$TMPD/e"
	ck
done

# -iobuf with -sorted and -nobuf together.
tc "-iobuf 64K -nobuf -sorted -wa -wb" "60000 records" "reference awk"
run -iobuf 64K -nobuf -sorted -a "$TMPD/biga" -b "$TMPD/bigb" -wa -wb
cp "$TMPD/bige" "$TMPD/e"
ck

# -iobuf has no effect on compressed input (help text): output is unchanged with 64K on gzipped input larger than the buffer.
tc "-iobuf 64K -wa -wb" "gzip 60000 records" "reference awk"
gzip -c "$TMPD/biga" > "$TMPD/biga.gz"
gzip -c "$TMPD/bigb" > "$TMPD/bigb.gz"
run -iobuf 64K -a "$TMPD/biga.gz" -b "$TMPD/bigb.gz" -wa -wb
cp "$TMPD/bige" "$TMPD/e"
ck

# A long line (70000 characters) with a small buffer.
tc "-iobuf 4096 -wa -wb" "bed+ with 70000-character column" "hand"
printf 'chr1\t10\t20\tname\t%s\n' "$LONG" > "$TMPD/a"
printf 'chr1\t15\t25\tn2\t%s\n' "$LONG" > "$TMPD/b"
run -iobuf 4096 -a "$TMPD/a" -b "$TMPD/b" -wa -wb
printf 'chr1\t10\t20\tname\t%s\tchr1\t15\t25\tn2\t%s\n' "$LONG" "$LONG" > "$TMPD/e"
ck

# -iobuf on stdin input.
tc "-iobuf 64K -a stdin" "60000 records on stdin" "reference awk"
RC=0; "$BT" intersect -iobuf 64K -a stdin -b "$TMPD/bigb" -wa -wb < "$TMPD/biga" > "$TMPD/o" 2>/dev/null || RC=$?
cp "$TMPD/bige" "$TMPD/e"
ck

# ---------------------------------------------------------------------
# -nobuf on small inputs: same bytes as buffered, derived by hand.
# ---------------------------------------------------------------------

printf 'chr1\t10\t20\tn1\t5\t+\nchr1\t30\t40\tn2\t7\t-\nchr1\t50\t60\tn3\t1\t+\n' > "$TMPD/a"
printf 'chr1\t15\t35\tx\t1\t+\nchr1\t100\t110\ty\t2\t-\n' > "$TMPD/b"

for m in "" "-wa" "-wb" "-wa -wb" "-wo" "-wao" "-loj" "-u" "-v" "-c"; do
	# -nobuf and then default buffering, mode '$m': both must equal the hand-written expectation.
	tc "-nobuf $m" "bed6" "hand"
	case "$m" in
	"") printf 'chr1\t15\t20\tn1\t5\t+\nchr1\t30\t35\tn2\t7\t-\n' ;;
	"-wa") printf 'chr1\t10\t20\tn1\t5\t+\nchr1\t30\t40\tn2\t7\t-\n' ;;
	"-wb") printf 'chr1\t15\t20\tn1\t5\t+\tchr1\t15\t35\tx\t1\t+\nchr1\t30\t35\tn2\t7\t-\tchr1\t15\t35\tx\t1\t+\n' ;;
	"-wa -wb") printf 'chr1\t10\t20\tn1\t5\t+\tchr1\t15\t35\tx\t1\t+\nchr1\t30\t40\tn2\t7\t-\tchr1\t15\t35\tx\t1\t+\n' ;;
	"-wo") printf 'chr1\t10\t20\tn1\t5\t+\tchr1\t15\t35\tx\t1\t+\t5\nchr1\t30\t40\tn2\t7\t-\tchr1\t15\t35\tx\t1\t+\t5\n' ;;
	"-wao") printf 'chr1\t10\t20\tn1\t5\t+\tchr1\t15\t35\tx\t1\t+\t5\nchr1\t30\t40\tn2\t7\t-\tchr1\t15\t35\tx\t1\t+\t5\nchr1\t50\t60\tn3\t1\t+\t.\t-1\t-1\t.\t-1\t.\t0\n' ;;
	"-loj") printf 'chr1\t10\t20\tn1\t5\t+\tchr1\t15\t35\tx\t1\t+\nchr1\t30\t40\tn2\t7\t-\tchr1\t15\t35\tx\t1\t+\nchr1\t50\t60\tn3\t1\t+\t.\t-1\t-1\t.\t-1\t.\n' ;;
	"-u") printf 'chr1\t10\t20\tn1\t5\t+\nchr1\t30\t40\tn2\t7\t-\n' ;;
	"-v") printf 'chr1\t50\t60\tn3\t1\t+\n' ;;
	"-c") printf 'chr1\t10\t20\tn1\t5\t+\t1\nchr1\t30\t40\tn2\t7\t-\t1\nchr1\t50\t60\tn3\t1\t+\t0\n' ;;
	esac > "$TMPD/e"
	run -nobuf -a "$TMPD/a" -b "$TMPD/b" $m
	ck
	tc "$m (buffered)" "bed6" "hand"
	run -a "$TMPD/a" -b "$TMPD/b" $m
	ck
done

# Buffered output of the same commands equals the same expected text (spot check on -wao and -c).
tc "-wao" "bed6" "hand"
printf 'chr1\t10\t20\tn1\t5\t+\tchr1\t15\t35\tx\t1\t+\t5\nchr1\t30\t40\tn2\t7\t-\tchr1\t15\t35\tx\t1\t+\t5\nchr1\t50\t60\tn3\t1\t+\t.\t-1\t-1\t.\t-1\t.\t0\n' > "$TMPD/e"
run -a "$TMPD/a" -b "$TMPD/b" -wao
ck

# -nobuf with -header.
tc "-nobuf -header" "bed6 with header" "hand"
{ printf '#h\n'; cat "$TMPD/a"; } > "$TMPD/ah"
run -nobuf -a "$TMPD/ah" -b "$TMPD/b" -header -u
printf '#h\nchr1\t10\t20\tn1\t5\t+\nchr1\t30\t40\tn2\t7\t-\n' > "$TMPD/e"
ck

# ---------------------------------------------------------------------
# Multiple databases. A has two chromosomes; the three databases are
# described in the comments on each expected block.
#   A (bed4): a1 chr1 10-30, a2 chr1 40-60, a3 chr1 70-90, a5 chr1 200-210, a4 chr2 10-30
#   DB1: chr1 20-25, chr1 50-55
#   DB2: chr1 25-45, chr1 100-110
#   DB3: chr1 28-29, chr1 80-85, chr2 5-15
# Hits per A, in database order then file order:
#   a1: DB1 20-25 (5 bases), DB2 25-45 (5), DB3 28-29 (1)
#   a2: DB1 50-55 (5), DB2 25-45 (5)
#   a3: DB3 80-85 (5)      a4: DB3 chr2 5-15 (5)      a5: none
# All coordinates are below 16 kb, so unsorted mode gives the same order.
# ---------------------------------------------------------------------

tb > "$TMPD/ma.bed" <<'EOF'
chr1|10|30|a1
chr1|40|60|a2
chr1|70|90|a3
chr1|200|210|a5
chr2|10|30|a4
EOF
tb > "$TMPD/m1.bed" <<'EOF'
chr1|20|25
chr1|50|55
EOF
tb > "$TMPD/m2.bed" <<'EOF'
chr1|25|45
chr1|100|110
EOF
tb > "$TMPD/m3.bed" <<'EOF'
chr1|28|29
chr1|80|85
chr2|5|15
EOF
MDB="$TMPD/m1.bed $TMPD/m2.bed $TMPD/m3.bed"

tb > "$TMPD/e.mwawb" <<'EOF'
chr1|10|30|a1|1|chr1|20|25
chr1|10|30|a1|2|chr1|25|45
chr1|10|30|a1|3|chr1|28|29
chr1|40|60|a2|1|chr1|50|55
chr1|40|60|a2|2|chr1|25|45
chr1|70|90|a3|3|chr1|80|85
chr2|10|30|a4|3|chr2|5|15
EOF
tb > "$TMPD/e.mwo" <<'EOF'
chr1|10|30|a1|1|chr1|20|25|5
chr1|10|30|a1|2|chr1|25|45|5
chr1|10|30|a1|3|chr1|28|29|1
chr1|40|60|a2|1|chr1|50|55|5
chr1|40|60|a2|2|chr1|25|45|5
chr1|70|90|a3|3|chr1|80|85|5
chr2|10|30|a4|3|chr2|5|15|5
EOF
tb > "$TMPD/e.mloj" <<'EOF'
chr1|10|30|a1|1|chr1|20|25
chr1|10|30|a1|2|chr1|25|45
chr1|10|30|a1|3|chr1|28|29
chr1|40|60|a2|1|chr1|50|55
chr1|40|60|a2|2|chr1|25|45
chr1|70|90|a3|3|chr1|80|85
chr1|200|210|a5|.|.|-1|-1
chr2|10|30|a4|3|chr2|5|15
EOF
tb > "$TMPD/e.mwao" <<'EOF'
chr1|10|30|a1|1|chr1|20|25|5
chr1|10|30|a1|2|chr1|25|45|5
chr1|10|30|a1|3|chr1|28|29|1
chr1|40|60|a2|1|chr1|50|55|5
chr1|40|60|a2|2|chr1|25|45|5
chr1|70|90|a3|3|chr1|80|85|5
chr1|200|210|a5|.|.|-1|-1|0
chr2|10|30|a4|3|chr2|5|15|5
EOF
tb > "$TMPD/e.mC" <<'EOF'
chr1|10|30|a1|1|1
chr1|10|30|a1|2|1
chr1|10|30|a1|3|1
chr1|40|60|a2|1|1
chr1|40|60|a2|2|1
chr1|40|60|a2|3|0
chr1|70|90|a3|1|0
chr1|70|90|a3|2|0
chr1|70|90|a3|3|1
chr1|200|210|a5|1|0
chr1|200|210|a5|2|0
chr1|200|210|a5|3|0
chr2|10|30|a4|1|0
chr2|10|30|a4|2|0
chr2|10|30|a4|3|1
EOF
tb > "$TMPD/e.mc" <<'EOF'
chr1|10|30|a1|3
chr1|40|60|a2|2
chr1|70|90|a3|1
chr1|200|210|a5|0
chr2|10|30|a4|1
EOF
tb > "$TMPD/e.mu" <<'EOF'
chr1|10|30|a1
chr1|40|60|a2
chr1|70|90|a3
chr2|10|30|a4
EOF
tb > "$TMPD/e.mv" <<'EOF'
chr1|200|210|a5
EOF
tb > "$TMPD/e.mwa" <<'EOF'
chr1|10|30|a1
chr1|10|30|a1
chr1|10|30|a1
chr1|40|60|a2
chr1|40|60|a2
chr1|70|90|a3
chr2|10|30|a4
EOF
tb > "$TMPD/e.mdef" <<'EOF'
chr1|20|25|a1
chr1|25|30|a1
chr1|28|29|a1
chr1|50|55|a2
chr1|40|45|a2
chr1|80|85|a3
chr2|10|15|a4
EOF
tb > "$TMPD/e.mwb" <<'EOF'
chr1|20|25|a1|1|chr1|20|25
chr1|25|30|a1|2|chr1|25|45
chr1|28|29|a1|3|chr1|28|29
chr1|50|55|a2|1|chr1|50|55
chr1|40|45|a2|2|chr1|25|45
chr1|80|85|a3|3|chr1|80|85
chr2|10|15|a4|3|chr2|5|15
EOF

for sorted in "" "-sorted"; do
	label=${sorted:-unsorted}

	# Three databases, -wa -wb: file id column 1, 2, 3 before each B record ($label).
	tc "-wa -wb -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -wa -wb $sorted
	cp "$TMPD/e.mwawb" "$TMPD/e"
	ck

	# Three databases, -wo: overlap length after the B record ($label).
	tc "-wo -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -wo $sorted
	cp "$TMPD/e.mwo" "$TMPD/e"
	ck

	# Three databases, -loj: a5 has a dummy '.' id and the BED3 NULL record ($label).
	tc "-loj -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -loj $sorted
	cp "$TMPD/e.mloj" "$TMPD/e"
	ck

	# Three databases, -wao ($label).
	tc "-wao -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -wao $sorted
	cp "$TMPD/e.mwao" "$TMPD/e"
	ck

	# Three databases, -C: one line per database per A with that database's hit count ($label).
	tc "-C -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -C $sorted
	cp "$TMPD/e.mC" "$TMPD/e"
	ck

	# Three databases, -c: one total count per A ($label).
	tc "-c -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -c $sorted
	cp "$TMPD/e.mc" "$TMPD/e"
	ck

	# Three databases, -u: A hit by any database, once ($label).
	tc "-u -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -u $sorted
	cp "$TMPD/e.mu" "$TMPD/e"
	ck

	# Three databases, -v: only a5 is hit by none ($label).
	tc "-v -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -v $sorted
	cp "$TMPD/e.mv" "$TMPD/e"
	ck

	# Three databases, -wa only: A repeated once per hit across all databases ($label).
	tc "-wa -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -wa $sorted
	cp "$TMPD/e.mwa" "$TMPD/e"
	ck

	# Three databases, no flags: A with the overlap coordinates for each hit, no id column ($label).
	tc "(none) -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB $sorted
	cp "$TMPD/e.mdef" "$TMPD/e"
	ck

	# Three databases, -wb only: overlap coordinates, then id and B ($label).
	tc "-wb -b 3 files $sorted" "3 dbs bed4/bed3" "hand"
	run -a "$TMPD/ma.bed" -b $MDB -wb $sorted
	cp "$TMPD/e.mwb" "$TMPD/e"
	ck
done

# -names replaces the id by the given tag, in database order.
tc "-wa -wb -names" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wa -wb -names one two three
awk -F'\t' -v OFS='\t' 'BEGIN { n[1] = "one"; n[2] = "two"; n[3] = "three" } { $5 = n[$5]; print }' "$TMPD/e.mwawb" > "$TMPD/e"
ck

# -names with -wo.
tc "-wo -names" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wo -names one two three
awk -F'\t' -v OFS='\t' 'BEGIN { n[1] = "one"; n[2] = "two"; n[3] = "three" } { $5 = n[$5]; print }' "$TMPD/e.mwo" > "$TMPD/e"
ck

# -names with -C: the tag replaces the database id.
tc "-C -names" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -C -names one two three
awk -F'\t' -v OFS='\t' 'BEGIN { n[1] = "one"; n[2] = "two"; n[3] = "three" } { $5 = n[$5]; print }' "$TMPD/e.mC" > "$TMPD/e"
ck

# -names with -loj: the dummy for a5 stays '.'.
tc "-loj -names" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -loj -names one two three
awk -F'\t' -v OFS='\t' 'BEGIN { n[1] = "one"; n[2] = "two"; n[3] = "three" } $5 in n { $5 = n[$5] } { print }' "$TMPD/e.mloj" > "$TMPD/e"
ck

# -names with -wa -wb in -sorted mode.
tc "-wa -wb -names -sorted" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wa -wb -names one two three -sorted
awk -F'\t' -v OFS='\t' 'BEGIN { n[1] = "one"; n[2] = "two"; n[3] = "three" } { $5 = n[$5]; print }' "$TMPD/e.mwawb" > "$TMPD/e"
ck

# -filenames prints the path as given on the command line.
tc "-wa -wb -filenames" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wa -wb -filenames
awk -F'\t' -v OFS='\t' -v d="$TMPD" '{ $5 = d "/m" $5 ".bed"; print }' "$TMPD/e.mwawb" > "$TMPD/e"
ck

# -filenames with -C.
tc "-C -filenames" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -C -filenames
awk -F'\t' -v OFS='\t' -v d="$TMPD" '{ $5 = d "/m" $5 ".bed"; print }' "$TMPD/e.mC" > "$TMPD/e"
ck

# -filenames with -wo, sorted.
tc "-wo -filenames -sorted" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wo -filenames -sorted
awk -F'\t' -v OFS='\t' -v d="$TMPD" '{ $5 = d "/m" $5 ".bed"; print }' "$TMPD/e.mwo" > "$TMPD/e"
ck

# -c with -names is rejected.
tc "-c -names" "3 dbs" "hand"
expect_error "-c cannot be used with -filenames or -names" "$BT" intersect -a "$TMPD/ma.bed" -b $MDB -c -names one two three

# Fewer names than databases is rejected.
tc "-names" "3 dbs" "hand"
expect_error "Number of database name tags given does not match" "$BT" intersect -a "$TMPD/ma.bed" -b $MDB -wa -wb -names one two

# -sortout orders the hits of each A by start then end: a1 keeps DB1, DB2, DB3 (starts 20, 25, 28)
# but a2 reports DB2 (25-45) before DB1 (50-55).
tc "-wa -wb -sortout" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wa -wb -sortout
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|chr1|20|25
chr1|10|30|a1|2|chr1|25|45
chr1|10|30|a1|3|chr1|28|29
chr1|40|60|a2|2|chr1|25|45
chr1|40|60|a2|1|chr1|50|55
chr1|70|90|a3|3|chr1|80|85
chr2|10|30|a4|3|chr2|5|15
EOF
ck

# -sortout with -sorted.
tc "-wa -wb -sortout -sorted" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wa -wb -sortout -sorted
ck

# -sortout with -names: the tags follow the reordered hits.
tc "-wa -wb -sortout -names" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wa -wb -sortout -names one two three
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|one|chr1|20|25
chr1|10|30|a1|two|chr1|25|45
chr1|10|30|a1|three|chr1|28|29
chr1|40|60|a2|two|chr1|25|45
chr1|40|60|a2|one|chr1|50|55
chr1|70|90|a3|three|chr1|80|85
chr2|10|30|a4|three|chr2|5|15
EOF
ck

# -sortout with -wo.
tc "-wo -sortout" "3 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB -wo -sortout
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|chr1|20|25|5
chr1|10|30|a1|2|chr1|25|45|5
chr1|10|30|a1|3|chr1|28|29|1
chr1|40|60|a2|2|chr1|25|45|5
chr1|40|60|a2|1|chr1|50|55|5
chr1|70|90|a3|3|chr1|80|85|5
chr2|10|30|a4|3|chr2|5|15|5
EOF
ck

# The same file given twice as two databases doubles every hit and numbers them 1 and 2.
tc "-wa -wb -b f f" "2 identical dbs" "hand"
run -a "$TMPD/ma.bed" -b "$TMPD/m1.bed" "$TMPD/m1.bed" -wa -wb
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|chr1|20|25
chr1|10|30|a1|2|chr1|20|25
chr1|40|60|a2|1|chr1|50|55
chr1|40|60|a2|2|chr1|50|55
EOF
ck

# ... and -C reports the same count for both.
tc "-C -b f f" "2 identical dbs" "hand"
run -a "$TMPD/ma.bed" -b "$TMPD/m1.bed" "$TMPD/m1.bed" -C
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|1
chr1|10|30|a1|2|1
chr1|40|60|a2|1|1
chr1|40|60|a2|2|1
chr1|70|90|a3|1|0
chr1|70|90|a3|2|0
chr1|200|210|a5|1|0
chr1|200|210|a5|2|0
chr2|10|30|a4|1|0
chr2|10|30|a4|2|0
EOF
ck

# Two databases, -C: ids 1 and 2.
tc "-C -b 2 files" "2 dbs" "hand"
run -a "$TMPD/ma.bed" -b "$TMPD/m1.bed" "$TMPD/m2.bed" -C
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|1
chr1|10|30|a1|2|1
chr1|40|60|a2|1|1
chr1|40|60|a2|2|1
chr1|70|90|a3|1|0
chr1|70|90|a3|2|0
chr1|200|210|a5|1|0
chr1|200|210|a5|2|0
chr2|10|30|a4|1|0
chr2|10|30|a4|2|0
EOF
ck

# One database, -C: no id column, a single count per A.
tc "-C -b 1 file" "1 db" "hand"
run -a "$TMPD/ma.bed" -b "$TMPD/m1.bed" -C
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1
chr1|40|60|a2|1
chr1|70|90|a3|0
chr1|200|210|a5|0
chr2|10|30|a4|0
EOF
ck

# Four databases: the fourth overlaps a3 and a5 only. -wo ids run 1 to 4.
tc "-wo -b 4 files" "4 dbs" "hand"
tb > "$TMPD/m4.bed" <<'EOF'
chr1|60|75
chr1|205|300
EOF
run -a "$TMPD/ma.bed" -b $MDB "$TMPD/m4.bed" -wo
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|chr1|20|25|5
chr1|10|30|a1|2|chr1|25|45|5
chr1|10|30|a1|3|chr1|28|29|1
chr1|40|60|a2|1|chr1|50|55|5
chr1|40|60|a2|2|chr1|25|45|5
chr1|70|90|a3|3|chr1|80|85|5
chr1|70|90|a3|4|chr1|60|75|5
chr1|200|210|a5|4|chr1|205|300|5
chr2|10|30|a4|3|chr2|5|15|5
EOF
ck

# Four databases, -C, with a4 untouched by DB1, DB2, DB4.
tc "-C -b 4 files" "4 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB "$TMPD/m4.bed" -C
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|1
chr1|10|30|a1|2|1
chr1|10|30|a1|3|1
chr1|10|30|a1|4|0
chr1|40|60|a2|1|1
chr1|40|60|a2|2|1
chr1|40|60|a2|3|0
chr1|40|60|a2|4|0
chr1|70|90|a3|1|0
chr1|70|90|a3|2|0
chr1|70|90|a3|3|1
chr1|70|90|a3|4|1
chr1|200|210|a5|1|0
chr1|200|210|a5|2|0
chr1|200|210|a5|3|0
chr1|200|210|a5|4|1
chr2|10|30|a4|1|0
chr2|10|30|a4|2|0
chr2|10|30|a4|3|1
chr2|10|30|a4|4|0
EOF
ck

# Four databases, -c.
tc "-c -b 4 files" "4 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB "$TMPD/m4.bed" -c
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|3
chr1|40|60|a2|2
chr1|70|90|a3|2
chr1|200|210|a5|1
chr2|10|30|a4|1
EOF
ck

# Four databases, -v: nothing is left over once DB4 covers a5.
tc "-v -b 4 files" "4 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB "$TMPD/m4.bed" -v
: > "$TMPD/e"
ck

# Four databases in -sorted mode, -wo.
tc "-wo -b 4 files -sorted" "4 dbs" "hand"
run -a "$TMPD/ma.bed" -b $MDB "$TMPD/m4.bed" -wo -sorted
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|chr1|20|25|5
chr1|10|30|a1|2|chr1|25|45|5
chr1|10|30|a1|3|chr1|28|29|1
chr1|40|60|a2|1|chr1|50|55|5
chr1|40|60|a2|2|chr1|25|45|5
chr1|70|90|a3|3|chr1|80|85|5
chr1|70|90|a3|4|chr1|60|75|5
chr1|200|210|a5|4|chr1|205|300|5
chr2|10|30|a4|3|chr2|5|15|5
EOF
ck

# Databases of different record types print each hit in its own format: DB1 bed6, DB2 bed4, DB3 gff.
tc "-wa -wb -b 3 mixed types" "3 dbs bed6/bed4/gff" "hand"
printf 'chr1\t10\t20\tp\t5\t+\n' > "$TMPD/x1.bed"
printf 'chr1\t15\t25\tq\n' > "$TMPD/x2.bed"
printf 'chr1\tsrc\tf\t18\t30\t.\t-\t.\tID=r;Note=a b\n' > "$TMPD/x3.gff"
printf 'chr1\t12\t22\tz\n' > "$TMPD/xa.bed"
run -a "$TMPD/xa.bed" -b "$TMPD/x1.bed" "$TMPD/x2.bed" "$TMPD/x3.gff" -wa -wb
tb > "$TMPD/e" <<'EOF'
chr1|12|22|z|1|chr1|10|20|p|5|+
chr1|12|22|z|2|chr1|15|25|q
chr1|12|22|z|3|chr1|src|f|18|30|.|-|.|ID=r;Note=a b
EOF
ck

# Same databases, -wo: overlaps 8, 7, 5 (GFF 18..30 is 0-based 17..30, overlap 17..22).
tc "-wo -b 3 mixed types" "3 dbs bed6/bed4/gff" "hand"
run -a "$TMPD/xa.bed" -b "$TMPD/x1.bed" "$TMPD/x2.bed" "$TMPD/x3.gff" -wo
tb > "$TMPD/e" <<'EOF'
chr1|12|22|z|1|chr1|10|20|p|5|+|8
chr1|12|22|z|2|chr1|15|25|q|7
chr1|12|22|z|3|chr1|src|f|18|30|.|-|.|ID=r;Note=a b|5
EOF
ck

# VCF A with BED databases: no id of the VCF, ids 1 and 2 for the databases, -wa -wb.
tc "-wa -wb -b 2 files" "vcf vs 2 dbs" "hand"
printf '##fileformat=VCFv4.2\nchr1\t150\tsnv1\tA\tG\t50\tPASS\t.\n' > "$TMPD/xv.vcf"
printf 'chr1\t100\t200\n' > "$TMPD/y1.bed"
printf 'chr1\t149\t150\tw\n' > "$TMPD/y2.bed"
run -a "$TMPD/xv.vcf" -b "$TMPD/y1.bed" "$TMPD/y2.bed" -wa -wb
tb > "$TMPD/e" <<'EOF'
chr1|150|snv1|A|G|50|PASS|.|1|chr1|100|200
chr1|150|snv1|A|G|50|PASS|.|2|chr1|149|150|w
EOF
ck

# Multiple databases with a header in A and -header.
tc "-header -wa -wb -b 3 files" "3 dbs with header in A" "hand"
{ printf '#hdr\n'; cat "$TMPD/ma.bed"; } > "$TMPD/mah.bed"
run -a "$TMPD/mah.bed" -b $MDB -header -wa -wb
{ printf '#hdr\n'; cat "$TMPD/e.mwawb"; } > "$TMPD/e"
ck

# Headers in every database are ignored.
tc "-wa -wb -b 3 files with headers" "3 dbs with headers" "hand"
for i in 1 2 3; do
	{ printf '#db%s\ntrack name=db%s\n' $i $i; cat "$TMPD/m$i.bed"; } > "$TMPD/mh$i.bed"
done
run -a "$TMPD/ma.bed" -b "$TMPD/mh1.bed" "$TMPD/mh2.bed" "$TMPD/mh3.bed" -wa -wb
cp "$TMPD/e.mwawb" "$TMPD/e"
ck

# An empty database file among three: ids stay 1, 2, 3 for the others.
tc "-wa -wb -b 3 files, one empty" "3 dbs, DB2 empty" "hand"
: > "$TMPD/m2empty.bed"
run -a "$TMPD/ma.bed" -b "$TMPD/m1.bed" "$TMPD/m2empty.bed" "$TMPD/m3.bed" -wa -wb
awk -F'\t' '$5 != 2' "$TMPD/e.mwawb" > "$TMPD/e"
ck

# A database with no records on A's chromosomes among three.
tc "-C -b 3 files, one on another chromosome" "3 dbs, DB2 other chromosome" "hand"
printf 'chr7\t1\t100\n' > "$TMPD/m2other.bed"
run -a "$TMPD/ma.bed" -b "$TMPD/m1.bed" "$TMPD/m2other.bed" "$TMPD/m3.bed" -C -nonamecheck
tb > "$TMPD/e" <<'EOF'
chr1|10|30|a1|1|1
chr1|10|30|a1|2|0
chr1|10|30|a1|3|1
chr1|40|60|a2|1|1
chr1|40|60|a2|2|0
chr1|40|60|a2|3|0
chr1|70|90|a3|1|0
chr1|70|90|a3|2|0
chr1|70|90|a3|3|1
chr1|200|210|a5|1|0
chr1|200|210|a5|2|0
chr1|200|210|a5|3|0
chr2|10|30|a4|1|0
chr2|10|30|a4|2|0
chr2|10|30|a4|3|1
EOF
ck

# -nobuf with three databases gives the -wa -wb text.
tc "-nobuf -wa -wb -b 3 files" "3 dbs" "hand"
run -nobuf -a "$TMPD/ma.bed" -b $MDB -wa -wb
cp "$TMPD/e.mwawb" "$TMPD/e"
ck

# -nobuf with three databases, -C.
tc "-nobuf -C -b 3 files" "3 dbs" "hand"
run -nobuf -a "$TMPD/ma.bed" -b $MDB -C
cp "$TMPD/e.mC" "$TMPD/e"
ck

# -nobuf, -sorted, -wao with three databases.
tc "-nobuf -sorted -wao -b 3 files" "3 dbs" "hand"
run -nobuf -sorted -a "$TMPD/ma.bed" -b $MDB -wao
cp "$TMPD/e.mwao" "$TMPD/e"
ck

# Databases of three kinds of compression with -C and -names.
tc "-C -names, gzip/bgzip/plain" "3 dbs compressed" "hand"
gzip -c "$TMPD/m1.bed" > "$TMPD/m1.bed.gz"
$htsutil bgzfcompress "$TMPD/m2.bed" "$TMPD/m2.bed.bgz"
run -a "$TMPD/ma.bed" -b "$TMPD/m1.bed.gz" "$TMPD/m2.bed.bgz" "$TMPD/m3.bed" -C -names p q r
awk -F'\t' -v OFS='\t' 'BEGIN { n[1] = "p"; n[2] = "q"; n[3] = "r" } { $5 = n[$5]; print }' "$TMPD/e.mC" > "$TMPD/e"
ck

# Gzipped A against three databases, -c.
tc "-c -b 3 files" "gzip A, 3 dbs" "hand"
gzip -c "$TMPD/ma.bed" > "$TMPD/ma.bed.gz"
run -a "$TMPD/ma.bed.gz" -b $MDB -c
cp "$TMPD/e.mc" "$TMPD/e"
ck

# ---------------------------------------------------------------------
# Multiple databases against a reference implementation (awk) on
# generated data: 600 A records and four databases of 400 records on
# chr1, chr10 and chr2, all sorted. Hits within an A are unordered in the
# comparison, so both sides go through sort.
# ---------------------------------------------------------------------

cat > "$TMPD/gen.awk" <<'EOF'
# Park-Miller generator, exact in double precision.
function rnd(n) { seed = (seed * 16807) % 2147483647; return seed % n }
BEGIN {
	seed = s0
	split("chr1 chr10 chr2", chrs, " ")
	for (i = 0; i < n; i++) {
		c = chrs[1 + rnd(3)]
		st = rnd(20000)
		len = 1 + rnd(400)
		printf "%s\t%d\t%d\n", c, st, st + len
	}
}
EOF
for k in 1 2 3 4; do
	awk -v s0=$((k * 7919)) -v n=400 -f "$TMPD/gen.awk" | sort -k1,1 -k2,2n -k3,3n > "$TMPD/r$k.bed"
done
awk -v s0=104729 -v n=600 -f "$TMPD/gen.awk" | sort -k1,1 -k2,2n -k3,3n > "$TMPD/ra.bed"
RDB="$TMPD/r1.bed $TMPD/r2.bed $TMPD/r3.bed $TMPD/r4.bed"

# ref <mode> <A> <db>...: reference implementation of the output of each mode.
cat > "$TMPD/ref.awk" <<'EOF'
FNR == NR { na++; A[na] = $0; ac[na] = $1; as[na] = $2; ae[na] = $3; next }
FNR == 1 { nd++ }
{ nb[nd]++; k = nb[nd]; B[nd, k] = $0; bc[nd, k] = $1; bs[nd, k] = $2; be[nd, k] = $3 }
END {
	for (i = 1; i <= na; i++) {
		total = 0
		for (d = 1; d <= nd; d++) {
			cnt[d] = 0
			for (k = 1; k <= nb[d]; k++) {
				if (bc[d, k] != ac[i] || !(as[i] < be[d, k] && bs[d, k] < ae[i])) continue
				cnt[d]++
				total++
				lo = as[i] > bs[d, k] ? as[i] : bs[d, k]
				hi = ae[i] < be[d, k] ? ae[i] : be[d, k]
				if (mode == "wawb") print A[i] "\t" d "\t" B[d, k]
				if (mode == "wo" || mode == "wao") print A[i] "\t" d "\t" B[d, k] "\t" (hi - lo)
				if (mode == "loj") print A[i] "\t" d "\t" B[d, k]
			}
		}
		if (total == 0) {
			if (mode == "loj") print A[i] "\t.\t.\t-1\t-1"
			if (mode == "wao") print A[i] "\t.\t.\t-1\t-1\t0"
			if (mode == "v") print A[i]
		} else if (mode == "u") print A[i]
		if (mode == "c") print A[i] "\t" total
		if (mode == "C") for (d = 1; d <= nd; d++) print A[i] "\t" d "\t" cnt[d]
	}
}
EOF

for m in wawb wo loj wao u v c C; do
	for sorted in "" "-sorted"; do
		# Four generated databases, mode -$m $sorted, compared as sorted lines with the awk reference.
		tc "-$m -b 4 files $sorted" "600 A records, 4 dbs of 400" "reference awk"
		case $m in
		wawb) flags="-wa -wb" ;;
		*) flags="-$m" ;;
		esac
		run -a "$TMPD/ra.bed" -b $RDB $flags $sorted
		sort "$TMPD/o" > "$TMPD/o.sorted"
		awk -v mode=$m -f "$TMPD/ref.awk" "$TMPD/ra.bed" $RDB | sort > "$TMPD/e"
		cksorted
	done
done

# The generated data compared with -nobuf and a 4K input buffer, -C -sorted.
tc "-nobuf -iobuf 4096 -C -sorted" "600 A records, 4 dbs of 400" "reference awk"
run -nobuf -iobuf 4096 -a "$TMPD/ra.bed" -b $RDB -C -sorted
sort "$TMPD/o" > "$TMPD/o.sorted"
awk -v mode=C -f "$TMPD/ref.awk" "$TMPD/ra.bed" $RDB | sort > "$TMPD/e"
cksorted

finish
