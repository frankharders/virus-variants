#!/usr/bin/env bash
# virus_ont_variants_coa.sh - Generiek virus Nanopore variant calling pipeline voor CoA rapportage
#
# Gebruik:
#   bash virus_ont_variants_coa.sh -r <referentie.fa> [opties]
#
# Verplicht:
#   -r FILE   Referentie FASTA
#
# Optioneel:
#   -g FILE   GFF3 annotatie (voor CSQ/gene/AA info in tabellen)
#   -i DIR    Map met Nanopore FASTQ reads (default: 01_basecalled)
#   -o DIR    Output map (default: results)
#   -t INT    Aantal threads (default: 8)
#   -d INT    Minimale read depth voor variant calling (default: 5)
#   -f FLOAT  Minimale allele frequentie voor variant calling (default: 0.05)
#   -m FLOAT  Minimale allele frequentie voor minority/CoA tabel (default: 0.05)
#   -q FLOAT  VarScan p-waarde drempel (default: 0.001)
#   -h        Toon deze helptext
#
# Output (per sample in <outdir>/tables/):
#   *.variants.tsv          Alle varianten >= -f, gesorteerd op positie
#   *.variants.extended.tsv Uitgebreide tabel met CSQ annotatie (alleen met -g)
#   *.minorities.tsv        Minority variants >= -m voor CoA rapportage
#   *.coa.tsv               CoA tabel met leesbare gen/AA namen
#   *.varscan_qc.tsv        AF distributie voor assay QC/afkapwaarde bepaling
#
# Variant calling: samtools mpileup | VarScan2 mpileup2snp
#   -B (geen BAQ correctie, aanbevolen voor ONT)
#   -Q 5 (lage base quality filter voor ONT)
#   --min-MQ 20 (mapping quality filter)
#   --strand-filter 0 (geen strand bias filter, niet van toepassing bij ONT)
#
# Vereiste tools: minimap2, samtools, bcftools, varscan, tabix, bgzip, awk, bc
#
# Auteur : WBVR Bioinformatics
# Versie : 3.0 (VarScan2 variant calling)

set -Eeuo pipefail

############################################
# Usage / argument parsing
############################################
usage(){
    awk 'NR==1{next} /^#!/{next} /^[^#]/{exit} {sub(/^# ?/,""); print}' "$0"
    exit "${1:-0}"
}

REF=""
GFF_IN=""
READS_DIR="01_basecalled"
OUT="results"
THREADS=8
MINREADS=5
MINAF=0.05
MINORITY_AF=0.05
VARSCAN_PVAL=0.001

while getopts ":r:g:i:o:t:d:f:m:q:h" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        g) GFF_IN="$OPTARG" ;;
        i) READS_DIR="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) MINREADS="$OPTARG" ;;
        f) MINAF="$OPTARG" ;;
        m) MINORITY_AF="$OPTARG" ;;
        q) VARSCAN_PVAL="$OPTARG" ;;
        h) usage 0 ;;
        :) echo "ERROR: Optie -$OPTARG vereist een argument." >&2; usage 1 ;;
        \?) echo "ERROR: Onbekende optie -$OPTARG." >&2; usage 1 ;;
    esac
done

# Validatie
errors=0
if [[ -z "$REF" ]]; then
    echo "ERROR: Referentie FASTA is verplicht (-r)." >&2; errors=1
fi
if [[ -n "$REF" && ! -f "$REF" ]]; then
    echo "ERROR: Referentie FASTA niet gevonden: $REF" >&2; errors=1
fi
if [[ -n "$GFF_IN" && ! -f "$GFF_IN" ]]; then
    echo "ERROR: GFF3 niet gevonden: $GFF_IN" >&2; errors=1
fi
if [[ ! -d "$READS_DIR" ]]; then
    echo "ERROR: Reads map niet gevonden: $READS_DIR" >&2; errors=1
fi
[[ $errors -gt 0 ]] && usage 1

# Veld-indexen voor Simple Table (bcftools csq output)
CSQ_FIELD_GENE=2
CSQ_FIELD_AA=6

############################################
# Helpers
############################################
ts(){ date +"%F %T"; }
say(){ echo "[$(ts)] $*" >&2; }
die(){ echo "[$(ts)] ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "$1 niet gevonden in PATH. Installeer $1 en probeer opnieuw."; }
hdr_id(){ head -1 "$1" | cut -d ' ' -f1 | sed 's/^>//'; }

sname(){
    bn=$(basename "$1")
    bn=${bn%%.fastq.gz}; bn=${bn%%.fq.gz}; bn=${bn%%.fastq}; bn=${bn%%.fq}
    echo "$bn"
}

# Extraheer FREQ uit VarScan FORMAT veld en converteer naar decimaal (0.1253)
varscan_af(){
    # $1 = FORMAT string (GT:GQ:SDP:DP:RD:AD:FREQ:...)
    # $2 = sample string
    awk -v fmt="$1" -v val="$2" 'BEGIN {
        split(fmt, f, ":")
        split(val, v, ":")
        for(i=1;i<=length(f);i++) if(f[i]=="FREQ") { gsub(/%/,"",v[i]); print v[i]/100; exit }
    }'
}

############################################
# Checks
############################################
say "STEP 0: Initialisatie en tool checks"
for t in minimap2 samtools bcftools varscan tabix bgzip awk sed grep sort bc; do need "$t"; done

say "  Referentie  : $REF"
say "  Reads map   : $READS_DIR"
say "  Output      : $OUT"
say "  Threads     : $THREADS"
say "  Min DP      : $MINREADS"
say "  Min AF      : $MINAF  (variant calling filter)"
say "  Minority AF : $MINORITY_AF  (CoA tabel drempel)"
say "  VarScan p   : $VARSCAN_PVAL"
[[ -n "$GFF_IN" ]] && say "  GFF3        : $GFF_IN" || say "  GFF3        : niet opgegeven (geen CSQ annotatie)"

############################################
# Init
############################################
mkdir -p "$OUT"
find "$OUT" -mindepth 1 -exec rm -rf -- {} + 2>/dev/null || true
mkdir -p "$OUT"/{alignments,variants,tables,logs,temp}
[[ -f "${REF}.fai" ]] || samtools faidx "$REF"

############################################
# GFF verwerken (stap 1)
############################################
CSQ_GFF=""
CSQ_OK=0

if [[ -n "$GFF_IN" ]]; then
    say "STEP 1: GFF-verwerking: normalisatie en synthese indien nodig"

    REFID="$(hdr_id "$REF")"
    GFF_NORM="$OUT/temp/input.relabel.gff3"

    say " - Normaliseren van GFF seqid naar FASTA ID: $REFID"
    AWK_SCRIPT='BEGIN{OFS="\t"} /^#/ {print; next} {$1=id; print}'
    if [[ "$GFF_IN" =~ \.gz$ ]]; then
        gzip -cd "$GFF_IN" | awk -v id="$REFID" "$AWK_SCRIPT" > "$GFF_NORM"
    else
        awk -v id="$REFID" "$AWK_SCRIPT" "$GFF_IN" > "$GFF_NORM"
    fi

    if ! grep -Eq $'\t(mRNA|transcript)\t' "$GFF_NORM"; then
        say " - Geen mRNA features gevonden: Synthetiseren van gene/mRNA/exon uit CDS"
        CSQ_GFF="$OUT/temp/input.csqready.gff3"
        awk -v OFS="\t" -v refid="$REFID" '
            BEGIN{ print "##gff-version 3" }
            function get_attr(h,key, i,n,a){
                n=split(h,a,";");
                for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]==key) return kv[2]}
                return ""
            }
            {
                c=$1; s=$2; f=$3; L=$4; R=$5; sc=$6; st=$7; ph=$8; at=$9;
                if (L ~ /^[0-9]+$/ && R ~ /^[0-9]+$/) {
                    if (f=="gene") {
                        gid=get_attr(at,"ID"); if(gid!=""){ g_strand[gid]=st; g_seen[gid]=1; g_name[gid]=gid }
                    }
                    else if (f=="CDS") {
                        g=get_attr(at,"gene"); if(g=="") g=get_attr(at,"Parent"); if(g=="") g=get_attr(at,"GeneID");
                        if(g==""){ g="gene" ++gcount }
                        cdsL[g][++idx[g]]=L; cdsR[g][idx[g]]=R; cdsSC[g][idx[g]]=(sc==""?"." : sc); cdsPH[g][idx[g]]=(ph==""?"." : ph);
                        cdsST[g]=st;
                        if(!(g in minL) || L<minL[g]) minL[g]=L;
                        if(!(g in maxR) || R>maxR[g]) maxR[g]=R;
                        g_name[g]=g; cds_seen[g]=1;
                    }
                    else { rest=rest $0 "\n" }
                }
            }
            END{
                for(g in cds_seen){
                    strand = (cdsST[g]!="") ? cdsST[g] : "+";
                    L=minL[g]; R=maxR[g];
                    gene_attrs = "ID=" g ";Name=" g ";gene=" g
                    print refid, "synthetic", "gene", L, R, ".", strand, ".", gene_attrs
                    mr = g ".t1"
                    mr_attrs = "ID=" mr ";Parent=" g ";Name=" mr ";gene=" g ";biotype=protein_coding"
                    print refid, "synthetic", "mRNA", L, R, ".", strand, ".", mr_attrs
                    nseg = asorti(cdsL[g], ord_idx)
                    for(i=1;i<=nseg;i++){
                        k = ord_idx[i];
                        eL = cdsL[g][k]; eR = cdsR[g][k];
                        exon_attrs = "ID=" mr ".exon" i ";Parent=" mr
                        print refid, "synthetic", "exon", eL, eR, ".", strand, ".", exon_attrs
                        cds_attrs = "Parent=" mr ";gene=" g
                        phase = cdsPH[g][k]; if(phase=="") phase="."
                        score = cdsSC[g][k]; if(score=="") score="."
                        print refid, "synthetic", "CDS", eL, eR, score, strand, phase, cds_attrs
                    }
                }
                if(rest!="") printf "%s", rest
            }
        ' "$GFF_NORM" > "$CSQ_GFF"

        if ! grep -Eq $'\tCDS\t' "$CSQ_GFF"; then
            say " - WAARSCHUWING: Synthese mislukt/levert geen CDS op. Annotatie wordt overgeslagen."
            CSQ_GFF=""
        else
            CSQ_OK=1
            say " - Synthese succesvol. Gebruik: $CSQ_GFF"
        fi
    else
        CSQ_GFF="$GFF_NORM"
        CSQ_OK=1
        say " - mRNA/transcript features gevonden. Gebruik: $CSQ_GFF"
    fi
else
    say "STEP 1: Geen GFF meegegeven → CSQ annotatie wordt overgeslagen (tabellen zonder Gene/AA info)"
fi

# Bouw gene-ID → leesbare naam lookup uit de GFF3 (voor CoA tabel)
GENE_LOOKUP_FILE="$OUT/temp/gene_lookup.tsv"
if [[ -n "$GFF_IN" ]]; then
    grep -v '^#' "$GFF_IN" | awk -F'\t' '$3=="gene"' | cut -f9 \
        | awk -F';' '{
            id=""; name="";
            for(i=1;i<=NF;i++){
                if($i~/^ID=/)   id=substr($i,4)
                if($i~/^Name=/) name=substr($i,6)
            }
            if(name=="") { name=id; sub(/^gene-/,"",name) }
            print id "\t" name " protein"
        }' > "$GENE_LOOKUP_FILE"
else
    touch "$GENE_LOOKUP_FILE"
fi

############################################
# Samples (stap 2 t/m 5)
############################################
shopt -s nullglob

READ_FILES=( "$READS_DIR"/*.fastq.gz "$READS_DIR"/*.fq.gz "$READS_DIR"/*.fastq "$READS_DIR"/*.fq )
TMP=(); for f in "${READ_FILES[@]}"; do [[ -f "$f" ]] && TMP+=("$f"); done; READ_FILES=( "${TMP[@]}" )
[[ ${#READ_FILES[@]} -gt 0 ]] || die "Geen FASTQs gevonden in $READS_DIR"

for READ_FILE in "${READ_FILES[@]}"; do

    SAMPLE="$(sname "$READ_FILE")"
    LOG="$OUT/logs/${SAMPLE}.log"

    say "STEP 2: [$SAMPLE] minimap2 aligneren → samtools sort/index"
    minimap2 -ax map-ont -t "$THREADS" "$REF" "$READ_FILE" 2> "$LOG" \
        | samtools sort -@ "$THREADS" -o "$OUT/alignments/${SAMPLE}.sorted.bam" - 2>> "$LOG"
    samtools index -@ "$THREADS" "$OUT/alignments/${SAMPLE}.sorted.bam" 2>> "$LOG"

    BAM_FILE="$OUT/alignments/${SAMPLE}.sorted.bam"
    RAW_VCF="$OUT/temp/${SAMPLE}.varscan.vcf"
    FINAL_VCF="$OUT/variants/${SAMPLE}.raw.vcf.gz"

    say "STEP 3: [$SAMPLE] VarScan2 variant calling"
    samtools mpileup -Q 5 -B -d 50000 --min-MQ 20 \
        -f "$REF" "$BAM_FILE" 2>> "$LOG" \
        | varscan mpileup2snp \
            --min-coverage "$MINREADS" \
            --min-var-freq "$MINAF" \
            --min-avg-qual 5 \
            --strand-filter 0 \
            --p-value "$VARSCAN_PVAL" \
            --output-vcf 1 \
        > "$RAW_VCF" 2>> "$LOG"

    # Comprimeer en indexeer
    bgzip -f "$RAW_VCF"
    tabix -f "${RAW_VCF}.gz"
    mv "${RAW_VCF}.gz" "$FINAL_VCF"
    mv "${RAW_VCF}.gz.tbi" "$FINAL_VCF.tbi"

    N_VARS=$(bcftools view -H "$FINAL_VCF" | wc -l)
    say "  -> $N_VARS varianten gevonden (>= ${MINAF})"

    CURRENT_SRC="$FINAL_VCF"
    VCF_CSQ_OK=0

    # --- STAP 4: CSQ annotatie ---
    if [[ "$CSQ_OK" -eq 1 ]]; then
        say "STEP 4: [$SAMPLE] bcftools csq annotatie"
        CSQ_VCF_RAW="$OUT/temp/${SAMPLE}.csq_raw.vcf.gz"
        CSQ_VCF="$OUT/variants/${SAMPLE}.csq.vcf.gz"
        if bcftools csq -f "$REF" -g "$CSQ_GFF" -p a -c CSQ -O z -o "$CSQ_VCF_RAW" "$CURRENT_SRC" >> "$LOG" 2>&1; then
            tabix -f "$CSQ_VCF_RAW"
            bcftools view "$CSQ_VCF_RAW" > "$OUT/temp/${SAMPLE}.csq_raw.vcf"
            say "  -> csq raw OK. Resolven van @pos MNV pointers..."
            python3 - "$OUT/temp/${SAMPLE}.csq_raw.vcf" "$CSQ_VCF" << 'PYEOF_INLINE'
import sys, re, gzip
in_file, out_file = sys.argv[1], sys.argv[2]
cache = {}
lines = open(in_file).readlines()
# Pass 1: build cache
for line in lines:
    if line.startswith("#"): continue
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 8: continue
    info = fields[7]
    if "CSQ=" in info and "CSQ=@" not in info:
        for part in info.split(";"):
            if part.startswith("CSQ="): cache[fields[1]] = part[4:]
# Pass 2: resolve and write
import subprocess
proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=open(out_file,"wb"))
for line in lines:
    if not line.startswith("#"):
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 8 and "CSQ=@" in fields[7]:
            new_info = []
            for part in fields[7].split(";"):
                if part.startswith("CSQ=@"):
                    ref_pos = part[5:]
                    part = "CSQ=" + cache.get(ref_pos, "@" + ref_pos)
                new_info.append(part)
            fields[7] = ";".join(new_info)
            line = "\t".join(fields) + "\n"
    proc.stdin.write(line.encode())
proc.stdin.close(); proc.wait()
PYEOF_INLINE
            tabix -f "$CSQ_VCF"
            rm -f "$OUT/temp/${SAMPLE}.csq_raw.vcf" "$CSQ_VCF_RAW" "$CSQ_VCF_RAW.tbi"
            VCF_CSQ_OK=1
            CURRENT_SRC="$CSQ_VCF"
            say "  -> csq OK (MNV pointers resolved). Gebruik $CSQ_VCF"
            say "  -> csq OK (MNV pointers resolved). Gebruik $CSQ_VCF"
        else
            say "  -> csq FAILED. Annotatie wordt overgeslagen (zie $LOG voor foutmelding)."
        fi
    fi

    say "STEP 5: [$SAMPLE] tabellen genereren"

    # VarScan FORMAT velden: GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR
    # FREQ is een percentage string (bijv "12.53%") — we zetten om naar decimaal in awk
    # DP = kolom 4 in FORMAT, FREQ = kolom 7

    # Gemeenschappelijke awk functie voor VarScan veld extractie
    VARSCAN_AWK_FIELDS='
        function get_freq(freq_str) {
            gsub(/%/, "", freq_str); return freq_str+0
        }
        function format_aa(aa_raw, consequence,    p, ref_aa, pos_aa, alt_aa, result) {
            # Geen AA change
            if (aa_raw == "") return ""
            # Missense formaten:
            # "1S>1P"   (pos+AA>pos+AA) → "S1P"
            # "15G>15D" (pos+AA>pos+AA) → "G15D"
            if (aa_raw ~ /^[0-9]+[A-Z*]>[0-9]+[A-Z*]$/) {
                split(aa_raw, p, ">")
                ref_aa = substr(p[1], length(p[1]), 1)
                pos_aa = substr(p[1], 1, length(p[1])-1)
                alt_aa = substr(p[2], length(p[2]), 1)
                result = ref_aa pos_aa alt_aa
            } else if (aa_raw ~ /^[A-Z*][0-9]+>[A-Z*][0-9]+$/) {
                # "L6373>Q6373" (AA+pos>AA+pos) → "L6373Q"
                split(aa_raw, p, ">")
                ref_aa = substr(p[1], 1, 1)
                pos_aa = substr(p[1], 2)
                alt_aa = substr(p[2], 1, 1)
                result = ref_aa pos_aa alt_aa
            } else {
                result = aa_raw
            }
            # Vervang * (stop codon) door X voor leesbaarheid in CoA rapportage
            gsub(/\*/, "X", result)
            return result
        }
    '

    # --- SIMPLE TABLE ---
    SIMPLE_TSV="$OUT/tables/${SAMPLE}.variants.tsv"
    say " - Simple table (gesorteerd op positie)"

    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        FMT_CSQ='%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t%INFO/CSQ\t1\n'
        FMT_NO_CSQ='%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n'
        {
            bcftools query -i 'INFO/CSQ!=""' -f "$FMT_CSQ" "$CURRENT_SRC"
            bcftools query -i 'INFO/CSQ=""'  -f "$FMT_NO_CSQ" "$CURRENT_SRC"
        } \
        | awk -v OFS="\t" -v G_IDX="$CSQ_FIELD_GENE" -v AA_IDX="$CSQ_FIELD_AA" \
              "$VARSCAN_AWK_FIELDS"'
            {
                pos=$1; ref=$2; alt=$3; dp=$4; freq_raw=$5; csq_str=$6; annotated=$7;
                freq = get_freq(freq_raw)
                gene=""; cDNA=""; aa_raw=""; consequence="";
                if (annotated == "1") {
                    split(csq_str, csq_vals, "|");
                    consequence=csq_vals[1]; gene=csq_vals[G_IDX]; aa_raw=csq_vals[AA_IDX];
                    cDNA = ref pos alt;
                    if (consequence ~ /synonymous/) aa_fmt = "Silent mutation"
                    else aa_fmt = format_aa(aa_raw, consequence)
                } else {
                    gene="Intergenic"; cDNA=ref pos alt; aa_fmt="";
                }
                if(length(ref)==1&&length(alt)==1){type="SNP";len=1}
                else if(length(ref)<length(alt)){type="INS";len=length(alt)-length(ref)}
                else if(length(ref)>length(alt)){type="DEL";len=length(ref)-length(alt)}
                else{type="VAR";len=length(alt)}
                vlabel=tolower(ref)""pos""tolower(alt);
                freq_str = sprintf("%.2f%%", freq)
                print type,vlabel,dp,len,freq_str,gene,cDNA,aa_fmt
            }' \
        | (
            echo -e "Variant Type\tVariant Position\tCoverage\tLength\tFrequency\tGene\tcDNA change\tAA change"
            sort -t $'\t' -k2,2V
        ) > "$SIMPLE_TSV"
    else
        bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\n' "$CURRENT_SRC" \
        | awk -v OFS="\t" "$VARSCAN_AWK_FIELDS"'
            {
                pos=$1; ref=$2; alt=$3; dp=$4; freq_raw=$5;
                freq = get_freq(freq_raw)
                if(length(ref)==1&&length(alt)==1){type="SNP";len=1}
                else if(length(ref)<length(alt)){type="INS";len=length(alt)-length(ref)}
                else if(length(ref)>length(alt)){type="DEL";len=length(ref)-length(alt)}
                else{type="VAR";len=length(alt)}
                vlabel=tolower(ref)""pos""tolower(alt);
                freq_str = sprintf("%.2f%%", freq)
                print type,vlabel,dp,len,freq_str,"Intergenic",ref pos alt,""
            }' \
        | (
            echo -e "Variant Type\tVariant Position\tCoverage\tLength\tFrequency\tGene\tcDNA change\tAA change"
            sort -t $'\t' -k2,2V
        ) > "$SIMPLE_TSV"
    fi

    # --- EXTENDED TABLE ---
    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        say " - Extended table (gesorteerd)"
        EXTENDED_TSV="$OUT/tables/${SAMPLE}.variants.extended.tsv"
        bcftools query -i 'INFO/CSQ!=""' \
            -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t%INFO/CSQ\t1\n' "$CURRENT_SRC" > "$OUT/temp/${SAMPLE}_csq.tmp"
        bcftools query -i 'INFO/CSQ=""' \
            -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC" >> "$OUT/temp/${SAMPLE}_csq.tmp"

        awk -v OFS="\t" "$VARSCAN_AWK_FIELDS"'
            {
                chrom=$1; pos=$2; ref=$3; alt=$4; dp=$5; freq_raw=$6; csq_str=$7; annotated=$8;
                if (annotated == "1") {
                    split(csq_str, f, "|");
                    consequence=f[1]; gene=f[2]; transcript=f[3]; biotype=f[4];
                    aa_raw=f[6]; split(f[7], dna_parts, ","); dna_change=dna_parts[1];
                    if (consequence ~ /synonymous/) aa_fmt = "Silent mutation"
                    else aa_fmt = format_aa(aa_raw, consequence)
                    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", \
                        chrom,pos,ref,alt,consequence,gene,transcript,biotype,dna_change,aa_fmt
                } else {
                    printf "%s\t%s\t%s\t%s\tIntergenic\tIntergenic\t\t\t\t\t\t%s\n", \
                        chrom,pos,ref,alt,ref pos alt
                }
            }' "$OUT/temp/${SAMPLE}_csq.tmp" \
        | (
            echo -e "CHROM\tPOS\tREF\tALT\tConsequence\tGene\tTranscript\tBiotype\tcDNA_change\tAA_change"
            sort -t $'\t' -k2,2n
        ) > "$EXTENDED_TSV"
        rm -f "$OUT/temp/${SAMPLE}_csq.tmp"
        say " - Extended table OK."
    fi

    # --- MINORITY VARIANTS TABLE (>= MINORITY_AF) ---
    MINORITY_PCT=$(echo "$MINORITY_AF * 100" | bc)
    say " - Minority variants table (AF >= ${MINORITY_PCT}%)"
    MINORITY_TSV="$OUT/tables/${SAMPLE}.minorities.tsv"

    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        {
            bcftools query -i 'INFO/CSQ!=""' \
                -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t%INFO/CSQ\t1\n' "$CURRENT_SRC"
            bcftools query -i 'INFO/CSQ=""' \
                -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
        }
    else
        bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
    fi \
    | awk -v OFS="\t" -v min_af="$MINORITY_AF" -v G_IDX="$CSQ_FIELD_GENE" -v AA_IDX="$CSQ_FIELD_AA" \
          "$VARSCAN_AWK_FIELDS"'
        {
            pos=$1; ref=$2; alt=$3; dp=$4; freq_raw=$5; csq_str=$6; annotated=$7;
            freq = get_freq(freq_raw)
            if (freq < min_af*100) next
            gene=""; aa_raw=""; consequence="";
            if (annotated == "1") {
                split(csq_str, csq_vals, "|");
                consequence = csq_vals[1]
                gene        = csq_vals[G_IDX];
                aa_raw      = csq_vals[AA_IDX];
            } else { gene = "Intergenic" }
            if(length(ref)==1 && length(alt)==1) type="SNP";
            else if(length(ref)<length(alt))     type="INS";
            else if(length(ref)>length(alt))     type="DEL";
            else                                  type="VAR";
            freq_str = sprintf("%.2f%%", freq)
            if (consequence ~ /synonymous/) aa_display = "Silent mutation"
            else aa_display = format_aa(aa_raw, consequence)
            print pos, type, ref, alt, freq_str, dp, gene, aa_display
        }' \
    | (
        echo -e "Position\tType\tRef\tAlt\tFrequency\tCoverage\tGene\tAA change"
        sort -t $'\t' -k1,1n
    ) > "$MINORITY_TSV"

    N_MIN=$(tail -n +2 "$MINORITY_TSV" | wc -l)
    say "  -> $N_MIN minority variant(s) >= ${MINORITY_PCT}% gevonden."

    # --- CoA TABEL ---
    say " - CoA tabel genereren"
    COA_TSV="$OUT/tables/${SAMPLE}.coa.tsv"

    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        {
            bcftools query -i 'INFO/CSQ!=""' \
                -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t%INFO/CSQ\t1\n' "$CURRENT_SRC"
            bcftools query -i 'INFO/CSQ=""' \
                -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
        }
    else
        bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\t\t0\n' "$CURRENT_SRC"
    fi \
    | awk -v OFS="\t" -v min_af="$MINORITY_AF" \
          -v G_IDX="$CSQ_FIELD_GENE" -v AA_IDX="$CSQ_FIELD_AA" \
          -v lookup_file="$GENE_LOOKUP_FILE" \
          "$VARSCAN_AWK_FIELDS"'
        BEGIN {
            while ((getline line < lookup_file) > 0) {
                split(line, kv, "\t")
                if (kv[1] != "") gene_name[kv[1]] = kv[2]
            }
            close(lookup_file)
        }
        {
            pos=$1; ref=$2; alt=$3; dp=$4; freq_raw=$5; csq_str=$6; annotated=$7;
            freq = get_freq(freq_raw)
            if (freq < min_af*100) next
            gene_id=""; aa_raw=""; consequence="";
            if (annotated == "1") {
                split(csq_str, csq_vals, "|");
                gene_id     = csq_vals[G_IDX];
                aa_raw      = csq_vals[AA_IDX];
                consequence = csq_vals[1];
            }
            # Gene display naam
            if (gene_id == "" || gene_id == "Intergenic") {
                gene_display = "Intergenic"
            } else if (gene_id in gene_name) {
                gene_display = gene_name[gene_id]
            } else {
                gene_display = gene_id; sub(/^gene-/, "", gene_display)
                gene_display = gene_display " protein"
            }
            # AA mutatie
            if (gene_id == "" || gene_id == "Intergenic") {
                aa_display = "Untranslated"
            } else if (consequence ~ /synonymous/) {
                aa_display = "Silent mutation"
            } else {
                aa_display = format_aa(aa_raw, consequence)
            }
            # Variant type en label
            if(length(ref)==1 && length(alt)==1) type="SNP";
            else if(length(ref)<length(alt))     type="INS";
            else if(length(ref)>length(alt))     type="DEL";
            else                                  type="VAR";
            vlabel = toupper(ref) pos toupper(alt)
            if(type=="SNP")      len=1
            else if(type=="INS") len=length(alt)-length(ref)
            else if(type=="DEL") len=length(ref)-length(alt)
            else                 len=length(alt)
            freq_str = sprintf("%.2f%%", freq)
            print type, vlabel, dp, len, freq_str, gene_display, aa_display
        }' \
    | (
        echo -e "Variant Type\tVariant Position and Identified Alternative Base\tCoverage\tLength of Variant\tFrequency of Variant\tGene (Region)\tAmino Acid Mutation"
        sort -t $'\t' -k2,2V
    ) > "$COA_TSV"

    N_COA=$(tail -n +2 "$COA_TSV" | wc -l)
    say "  -> CoA tabel: $N_COA variant(s) >= ${MINORITY_PCT}%."

    # --- QC DISTRIBUTIE TABEL ---
    say " - QC distributie tabel (AF verdeling voor assay kwaliteitsbepaling)"
    QC_TSV="$OUT/tables/${SAMPLE}.varscan_qc.tsv"
    QC_RAW_VCF="$OUT/temp/${SAMPLE}.varscan_qc_raw.vcf"
    samtools mpileup -Q 5 -B -d 50000 --min-MQ 20 \
        -f "$REF" "$BAM_FILE" 2>> "$LOG" \
        | varscan mpileup2snp \
            --min-coverage "$MINREADS" \
            --min-var-freq 0.02 \
            --min-avg-qual 5 \
            --strand-filter 0 \
            --p-value "$VARSCAN_PVAL" \
            --output-vcf 1 \
        > "$QC_RAW_VCF" 2>> "$LOG"


    bcftools query -f '%POS\t%REF\t%ALT\t[%DP]\t[%FREQ]\n' "$QC_RAW_VCF" \
    | awk "$VARSCAN_AWK_FIELDS"'
        {
            freq = get_freq($5)
            if      (freq <  3) bin="02-03%"
            else if (freq <  5) bin="03-05%"
            else if (freq < 10) bin="05-10%"
            else if (freq < 20) bin="10-20%"
            else if (freq < 50) bin="20-50%"
            else                bin=">50%  "
            count[bin]++
            total++
        }
        END {
            print "AF bin\tAantal varianten\tPercentage van totaal"
            bins[1]="02-03%"; bins[2]="03-05%"; bins[3]="05-10%"
            bins[4]="10-20%"; bins[5]="20-50%"; bins[6]=">50%  "
            for(i=1;i<=6;i++){
                b=bins[i]; n=(b in count)?count[b]:0
                printf "%s\t%d\t%.1f%%\n", b, n, (total>0?n/total*100:0)
            }
            printf "Totaal\t%d\t100.0%%\n", total
        }' > "$QC_TSV"

    rm -f "$QC_RAW_VCF"
    rm -f "$QC_RAW_VCF"
    say "  -> QC tabel geschreven. Zie $QC_TSV"
    say "DONE: [$SAMPLE]"
done

say "ALL DONE. Resultaten in $OUT. Gebruik IGV met $OUT/alignments/*.sorted.bam en $REF"