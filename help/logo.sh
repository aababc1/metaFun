
#!/bin/bash

# Directory to save logos
LOGO_DIR="../help"

# Generate logo for each module
figlet -f slant "INTERACTIVE 
MODULE for
COMPARATIVE_
ANNOTATION" > "${LOGO_DIR}/slant_logo_imca.txt"

figlet -f slant "INTERACTIVE 
MODULE for
WMS_
TAXONOMY" > "${LOGO_DIR}/slant_logo_imwt.txt"

figlet -f slant "RAWREAD_
QC" > "${LOGO_DIR}/slant_logo_rawread_qc.txt"

figlet -f slant "ASSEMBLY_
BINNING" > "${LOGO_DIR}/slant_logo_assembly_binning.txt"

figlet -f slant "BIN_
ASSESSMENT" > "${LOGO_DIR}/slant_logo_bin_assessment.txt"

figlet -f slant "WMS_
FUNCTION" > "${LOGO_DIR}/slant_logo_wms_function.txt"

figlet -f slant "WMS_
TAXONOMY" > "${LOGO_DIR}/slant_logo_wms_taxonomy.txt"

figlet -f slant "COMPARATIVE_
ANNOTATION" > "${LOGO_DIR}/slant_logo_comparative_annotation.txt"

figlet -f slant "GENOME
SELECTOR" > "${LOGO_DIR}/slant_logo_genome_selector.txt"
