# >>>>> YES
## TAXONOMY TRAINING
## Only needed ONCE per database

# UNITE
# ITS database: https://unite.ut.ee/repository.php
# Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2020): UNITE QIIME release for Fungi. Version 04.02.2020. UNITE Community. https://doi.org/10.15156/BIO/786385
# Extract from the reference database the amplicon that was amplified in your study
# FUNGAL note: In our experience, fungal ITS classifiers trained on the UNITE reference database do NOT benefit from extracting/trimming reads to primer sites. We recommend training UNITE classifiers on the full reference sequences. Furthermore, we recommend the “developer” sequences (located within the QIIME-compatible release download) because the standard versions of the sequences have already been trimmed to the ITS region (excluding portions of flanking rRNA genes that may be present in amplicons generated with standard ITS primers).
# HOWEVER: When I try to use the "developer" version of the database, the training fails with the message "Invalid character in sequence". You need to clean the sequences as instructed below.
# HOWEVER: In some other texts the authors claim that extraction with primers sometimes improves the classification.
# NOTE: We use the larger database, which "Includes global and 97% singletons.", because it outputs a better classification (less unassigned sequences)

DBPATH=/home/nezapa/qiime-thesis/database/200204UNITE/developer

# # if you use the "developer" version of UNITE, remove the non-IUPAC sequences from the dataset
# awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' $DBPATH/sh_refs_qiime_ver8_dynamic_s_04.02.2020_dev.fasta | sed -e '/^>/!s/\(.*\)/\U\1/;s/[[:blank:]]*$//' > $DBPATH/sh_refs_qiime_ver8_dynamic_s_04.02.2020.fasta
# ln -s $DBPATH/sh_taxonomy_qiime_ver8_dynamic_s_04.02.2020_dev.txt $DBPATH/sh_taxonomy_qiime_ver8_dynamic_s_04.02.2020.txt

# Import the reference database if it is not yet in QIIME format
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path $DBPATH/sh_refs_qiime_ver8_dynamic_s_04.02.2020.fasta \
  --output-path $DBPATH/dynamic_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path $DBPATH/sh_taxonomy_qiime_ver8_dynamic_s_04.02.2020.txt \
  --output-path $DBPATH/dynamic_taxonomy.qza

# remove the sequences with poor taxonomic classification (sequences named as unclassified - if these are the most similar to your sequence, they will end up as the classification of the ASV even if there is another very similar database sequence with a much better taxonomic classification; the taxonomy file does not need to be filtered; use a comma to delimit different filter expressions
qiime taxa filter-seqs \
	--i-sequences $DBPATH/dynamic_otus.qza \
	--i-taxonomy $DBPATH/dynamic_taxonomy.qza \
	--p-exclude g__unidentified \
	--o-filtered-sequences $DBPATH/dynamic_otus.filtered.qza

# filter taxonomy (even though it is supposedly not necessary)
grep -v "g__unidentified" $DBPATH/sh_taxonomy_qiime_ver8_dynamic_s_04.02.2020.txt > $DBPATH/dynamic_taxonomy.filtered.txt
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path $DBPATH/dynamic_taxonomy.filtered.txt \
  --output-path $DBPATH/dynamic_taxonomy.filtered.qza

# not neccessarily good for ITS, see recommendation above
qiime feature-classifier extract-reads \
  --i-sequences $DBPATH/dynamic_otus.filtered.qza \
  --p-f-primer ACTTTYRRCAAYGGATCWCT \
  --p-r-primer GCCTCCGCTTATTGATATGCTTAART \
  --o-reads $DBPATH/ref-seqs_58SFun_ITS4Fun.qza

# Train the classifier with full length sequences
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $DBPATH/dynamic_otus.filtered.qza \
  --i-reference-taxonomy $DBPATH/dynamic_taxonomy.filtered.qza \
  --o-classifier $DBPATH/whole.length.classifier.qza

# Train the classifier with primer-extracted sequences
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $DBPATH/ref-seqs_58SFun_ITS4Fun.qza \
  --i-reference-taxonomy $DBPATH/dynamic_taxonomy.qza \
  --o-classifier $DBPATH/classifier_58SFun_ITS4Fun.qza


# SILVA
# Latest SILVA database https://www.arb-silva.de/download/archive/qiime
# SILVA does not officially support QIIME at this moment. The QIIME team may provide updated files on their project site: http://qiime.org/home_static/dataFiles.html
# ALSO check this for possibly most recent database: https://docs.qiime2.org/2020.11/data-resources/
# I use the latter, full-length
# Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Opens external link in new windowNucl. Acids Res. 41 (D1): D590-D596.
# Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Opens external link in new windowNucl. Acids Res. 42:D643-D648
# wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip

DBPATH=/home/nezapa/qiime-thesis/database/200827SILVA138

# Train the classifier with complete sequences
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $DBPATH/silva-138-99-seqs.qza \
  --i-reference-taxonomy $DBPATH/silva-138-99-tax.qza \
  --o-classifier $DBPATH/classifier.all.qza


# extract the relevant part of the sequences - BACTERIA
qiime feature-classifier extract-reads \
  --i-sequences $DBPATH/silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --o-reads $DBPATH/ref-seqs_341F_802R.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $DBPATH/ref-seqs_341F_802R.qza \
  --i-reference-taxonomy $DBPATH/silva-138-99-tax.qza \
  --o-classifier $DBPATH/classifier-B-341F_802R.qza


# extract the relevant part of the sequences - ARCHAEA
qiime feature-classifier extract-reads \
  --i-sequences $DBPATH/silva-138-99-seqs.qza \
  --p-f-primer TCCGGTTGATCCYGCBRG \
  --p-r-primer GCTACGRRYGYTTTARRC \
  --o-reads $DBPATH/ref-seqs_SSU1ArF_SSU520R.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $DBPATH/ref-seqs_SSU1ArF_SSU520R.qza \
  --i-reference-taxonomy $DBPATH/silva-138-99-tax.qza \
  --o-classifier $DBPATH/classifier-A-SSU1ArF_SSU520R.qza


# extract the relevant part of the sequences - EUKARYA
qiime feature-classifier extract-reads \
  --i-sequences $DBPATH/silva-138-99-seqs.qza \
  --p-f-primer TTAAARVGYTCGTAGTYG \
  --p-r-primer CCGTCAATTHCTTYAART \
  --o-reads $DBPATH/ref-seqs_18S616F_18S1132R.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $DBPATH/ref-seqs_18S616F_18S1132R.qza \
  --i-reference-taxonomy $DBPATH/silva-138-99-tax.qza \
  --o-classifier $DBPATH/classifier-E-18S616F_18S1132R.qza

# ### END TAXONOMY TRAINING