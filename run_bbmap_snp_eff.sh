#!/usr/bin/env bash
# Nieuwe versie V7: Inclusief alle varianten (met en zonder CSQ-annotatie), gesorteerd op Variant Position
set -Eeuo pipefail

############################################
# Config
############################################
REF="references/MN908947.3.fa"
READS_DIR="RAWREADS"
OUT="results"
THREADS="${THREADS:-8}"
GFF_IN="${1:-}"

# BBTools CallVariants parameters
PLOIDY=1
MINREADS=5
MINAF=0.02
CALL_OPTS="ploidy=${PLOIDY} minreads=${MINREADS} minallelefraction=${MINAF}"

# Dit zijn de veld-indexen voor de Simple Table, gebaseerd op bcftools csq output
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
sname(){ bn=$(basename "$1"); bn=${bn%%.fastq.gz}; bn=${bn%%.fq.gz}; bn=${bn%%.fastq}; bn=${bn%%.fq}; echo "${bn/_R1/}"; }

############################################
# Checks
############################################
say "STEP 0: Initialisatie en tool checks"
for t in bbmap.sh callvariants.sh samtools bcftools tabix bgzip awk sed grep sort; do need "$t"; done
[[ -f "$REF" ]] || die "Referentie FASTA ontbreekt: $REF"

# ... (De rest van de initialisatie en GFF-verwerking (Stap 1) blijft ongewijzigd) ...

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
    [[ -f "$GFF_IN" ]] || die "GFF niet gevonden: $GFF_IN"
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
        # Awk-logica voor synthese
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


############################################
# Samples (stap 2 t/m 5)
############################################
shopt -s nullglob
R1S=( "$READS_DIR"/*_R1*.fastq.gz "$READS_DIR"/*_R1*.fq.gz "$READS_DIR"/*_R1*.fastq "$READS_DIR"/*_R1*.fq )
TMP=(); for f in "${R1S[@]}"; do [[ -f "$f" ]] && TMP+=("$f"); done; R1S=( "${TMP[@]}" )
[[ ${#R1S[@]} -gt 0 ]] || die "Geen R1 FASTQs in $READS_DIR"

for R1 in "${R1S[@]}"; do
    R2="${R1/_R1/_R2}"
    [[ -f "$R2" ]] || { say "WARN: [$R1] mist bijbehorende R2 → overslaan"; continue; }
    SAMPLE="$(sname "$R1")"
    LOG="$OUT/logs/${SAMPLE}.log"

    say "STEP 2: [$SAMPLE] BBMap aligneren → samtools sort/index"
    bbmap.sh ref="$REF" in="$R1" in2="$R2" outm=stdout.sam ambiguous=random maxindel=2000 local=t nodisk=t threads="$THREADS" 2> "$LOG" \
        | samtools sort -@ "$THREADS" -o "$OUT/alignments/${SAMPLE}.sorted.bam" - 2>> "$LOG"
    samtools index -@ "$THREADS" "$OUT/alignments/${SAMPLE}.sorted.bam" 2>> "$LOG"

    say "STEP 3: [$SAMPLE] callvariants.sh"
    RAW_VCF="$OUT/variants/${SAMPLE}.raw.vcf"
    callvariants.sh in="$OUT/alignments/${SAMPLE}.sorted.bam" ref="$REF" vcf="$RAW_VCF" $CALL_OPTS >> "$LOG" 2>&1
    bgzip -f "$RAW_VCF"; tabix -f "${RAW_VCF}.gz"

    CURRENT_SRC="${RAW_VCF}.gz"
    VCF_CSQ_OK=0 
    
    if [[ "$CSQ_OK" -eq 1 ]]; then
        say "STEP 4: [$SAMPLE] bcftools csq annotatie"
        CSQ_VCF="$OUT/variants/${SAMPLE}.csq.vcf.gz"
        if bcftools csq -f "$REF" -g "$CSQ_GFF" -p a -c CSQ -O z -o "$CSQ_VCF" "${RAW_VCF}.gz" >> "$LOG" 2>&1; then
            tabix -f "$CSQ_VCF"
            CURRENT_SRC="$CSQ_VCF"
            VCF_CSQ_OK=1
            say "  -> csq OK. Gebruik $CSQ_VCF"
        else
            say "  -> csq FAILED. Annotatie wordt overgeslagen (zie $LOG voor foutmelding)."
        fi
    fi

    say "STEP 5: [$SAMPLE] tabellen genereren"
    
    # --- SIMPLE TABLE: FIX VOOR INCLUSIE VAN NIET-GEANNOTTEERDE VARIANTEN EN SORTERING ---
    SIMPLE_TSV="$OUT/tables/${SAMPLE}.variants.tsv"
    
    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        say " - Simple table (met directe CSQ-parsing, Intergenic call en Sortering op Positie)"
        
        # 1. QUERY VOOR GEANNOTTEERDE VARIANTEN:
        FMT_CSQ="%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\t%INFO/CSQ\t1\n"
        
        # 2. QUERY VOOR NIET-GEANNOTTEERDE VARIANTEN:
        FMT_NO_CSQ="%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\t\t\t0\n"

        # Combineer beide query's en verwerk ze met AWK en SORT
        {
            bcftools query -i 'INFO/CSQ!=""' -f "$FMT_CSQ" "$CURRENT_SRC"
            bcftools query -e 'INFO/CSQ!=""' -f "$FMT_NO_CSQ" "$CURRENT_SRC"
        } \
        | awk -v OFS="\t" -v G_IDX="$CSQ_FIELD_GENE" -v AA_IDX="$CSQ_FIELD_AA" '
            {
                pos=$1; ref=$2; alt=$3; dp=$4; af=$5; csq_str=$6; annotated=$7;
                
                gene=""; cDNA=""; aa="";
                
                if (annotated == "1") {
                    split(csq_str, csq_vals, "|");
                    gene=csq_vals[G_IDX];
                    aa=csq_vals[AA_IDX];
                    
                    cDNA = "c." ref pos alt;
                    if (tolower(gene) !~ /orf/ && tolower(gene) !~ /spike/) {
                        cDNA = ref pos alt; 
                    }
                } else {
                    gene="Intergenic";
                    cDNA=ref pos alt;
                    aa="";
                }

                if(length(ref)==1&&length(alt)==1){type="SNP";len=1}
                else if(length(ref)<length(alt)){type="INS";len=length(alt)-length(ref)}
                else if(length(ref)>length(alt)){type="DEL";len=length(ref)-length(alt)}
                else{type="VAR";len=length(alt)}
                vlabel=tolower(ref)""pos""tolower(alt);
                
                split(af, af_vals, ","); 
                freq=(af_vals[1]==""?"":sprintf("%.2f%%",100*af_vals[1]));
                
                print type,vlabel,dp,len,freq,gene,cDNA,aa
            }' \
        | (
            # Voeg de header toe en sorteer de data
            echo -e "Variant Type\tVariant Position\tCoverage\tLength\tFrequency\tGene\tcDNA change\tAA change"
            # Sorteer op de 2e kolom (Variant Position) als numerieke waarde (V)
            sort -t $'\t' -k2,2V
        ) > "$SIMPLE_TSV"
            
    else
        say " - Simple table (fallback zonder annotaties, gesorteerd)"
        # Fallback zonder CSQ/BCSQ (nu ook gesorteerd)
        bcftools query -f "%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n" "$CURRENT_SRC" \
        | awk -v OFS="\t" '
            {
                pos=$1; ref=$2; alt=$3; dp=$4; af=$5;
                if(length(ref)==1&&length(alt)==1){type="SNP";len=1}
                else if(length(ref)<length(alt)){type="INS";len=length(alt)-length(ref)}
                else if(length(ref)>length(alt)){type="DEL";len=length(ref)-length(alt)}
                else{type="VAR";len=length(alt)}
                vlabel=tolower(ref)""pos""tolower(alt);
                
                gene="Intergenic";
                cDNA=ref pos alt;
                aa="";
                
                split(af, af_vals, ","); 
                freq=(af_vals[1]==""?"":sprintf("%.2f%%",100*af_vals[1]));
                
                print type,vlabel,dp,len,freq,gene,cDNA,aa
            }' \
        | (
            echo -e "Variant Type\tVariant Position\tCoverage\tLength\tFrequency\tGene\tcDNA change\tAA change"
            sort -t $'\t' -k2,2V
        ) > "$SIMPLE_TSV"
    fi

    # --- EXTENDED TABLE (Nu ook gesorteerd) ---
    if [[ "$VCF_CSQ_OK" -eq 1 ]]; then
        say " - Extended table (gesorteerd)"
        
        EXTENDED_TSV="$OUT/tables/${SAMPLE}.variants.extended.tsv"

        HDR_CSQ='CHROM\tPOS\tREF\tALT\tConsequence\tGene\tTranscript\tBiotype\tProtein_position\tAmino_acids\tCodons\tcDNA_change'
        FMT_CSQ_EXT='%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ[1]\t%INFO/CSQ[2]\t%INFO/CSQ[3]\t%INFO/CSQ[4]\t%INFO/CSQ[7]\t%INFO/CSQ[6]\t%INFO/CSQ[5]\t%POS\t%REF\t%ALT\t1\n' # '1' = annotated
        FMT_NO_CSQ_EXT='%CHROM\t%POS\t%REF\t%ALT\tIntergenic\tIntergenic\t\t\t\t\t\t%POS\t%REF\t%ALT\t0\n' # '0' = not annotated

        {
            bcftools query -i 'INFO/CSQ!=""' -f "$FMT_CSQ_EXT" "$CURRENT_SRC"
            bcftools query -e 'INFO/CSQ!=""' -f "$FMT_NO_CSQ_EXT" "$CURRENT_SRC"
        } \
        | awk -v OFS="\t" '
            {
                # $1-$4 (basis), $5-$11 (CSQ velden), $12-$14 (POS, REF, ALT voor cDNA), $15 (annotated flag)
                
                # Genereer de cDNA change
                cDNA_notatie = $13 $12 $14; # REF POS ALT
                
                # Verwijder de laatste 4 kolommen ($12 t/m $15) en voeg de cDNA notatie toe
                printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, cDNA_notatie
            }' \
            | (
                echo -e "$HDR_CSQ"
                # Sorteer op de 2e kolom (POS) als numerieke waarde (n)
                sort -t $'\t' -k2,2n
            ) > "$OUT/tables/${SAMPLE}.variants.extended.tsv.tmp" \
            && mv "$OUT/tables/${SAMPLE}.variants.extended.tsv.tmp" "$EXTENDED_TSV"
        
        say " - Extended table OK."
    fi

    say "DONE: [$SAMPLE]"
done

say "ALL DONE. Resultaten in $OUT. Gebruik IGV met $OUT/alignments/*.sorted.bam en $REF"
