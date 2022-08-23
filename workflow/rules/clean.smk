rule cleanAfterAlignment:
    shell: """
find results/*/*/*/  -maxdepth 0 -not -name "bam" -exec rm -rf {} \;
rm -rf results/trackhub    
"""
