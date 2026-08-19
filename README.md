# mbgd-pipeline
MBGD pipeline for orhology analysis
### Initial setup for programs and environments
      ./setup.sh
### Database setup for orthology assignment 
      ./setup_db.sh

### Setup before running
      source ${MPL_TOPDIR}/etc/bashrc

### Execlute orthology assignment based on mmseqs profile search
      exec_assignment.sh query_genome.fa [profdb]

### Execlute orthology clustering 

All proteome seqeuces in FAST format (.faa) and annotations in GFF format (.gff) are assumed to be stored in ```in_data``` directory.   

      exec_clustering.sh [file_prefix]
